"""One test per rule of Plan.md "Rules enforced" / "Custom annotation sidecars" (offline)."""
import pytest
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SimpleLocation, WithinPosition

from conftest import ID, random_dna, tag
from gbk_creator import (AnnotationError, FeatureError, GbkCreatorError, GenBankBuilder, NamingError,
                         TaxonomyError)


def loc(start, end, strand=1):
    return SimpleLocation(start, end, strand=strand)


def cds(b, n, start, end, strand=1, contig=None, **kw):
    kw.setdefault("product", "hypothetical protein")
    kw.setdefault("translation", "MKV")
    b.add_gene(tag(n), contig or f"{ID}_scf1", loc(start, end, strand), "CDS", **kw)


def trna(b, n, start, end, strand=1, **kw):
    b.add_gene(tag(n), f"{ID}_scf1", loc(start, end, strand), "tRNA", product="tRNA-Gly", **kw)


# ----------------------------------------------------------------------------- naming
class TestIdentifierAndStrain:
    def test_identifier_max_50_chars(self, make_builder):
        make_builder("A" * 50, strain="A")
        with pytest.raises(NamingError, match="51 characters"):
            make_builder("A" * 51, strain="A")

    @pytest.mark.parametrize("bad", ["a b", "a:b", "a/b", "a\\b", "a\tb", "a|b", "ä"])
    def test_identifier_forbidden_characters(self, make_builder, bad):
        with pytest.raises(NamingError, match="forbidden characters"):
            make_builder(bad, strain=None)

    def test_identifier_empty(self, make_builder):
        with pytest.raises(NamingError):
            make_builder("", strain=None)

    def test_strain_defaults_to_identifier(self, make_builder):
        assert make_builder(strain=None).strain == ID

    def test_strain_may_equal_identifier(self, make_builder):
        assert make_builder(strain=ID).strain == ID

    def test_identifier_must_start_with_strain(self, make_builder):
        with pytest.raises(NamingError, match="must start with the strain"):
            make_builder(strain="OTHER")
        with pytest.raises(NamingError, match="must start with the strain"):
            make_builder(strain="496-i1")  # a suffix, not a prefix

    def test_strain_max_40_chars(self, make_builder):
        make_builder("S" * 45, strain="S" * 40)
        with pytest.raises(NamingError, match="max 40"):
            make_builder("S" * 45, strain="S" * 41)

    @pytest.mark.parametrize("taxid", [0, -5, "1308", 1.5, True])
    def test_taxid_must_be_positive_int(self, make_builder, taxid):
        with pytest.raises(GbkCreatorError, match="taxid"):
            make_builder(taxid=taxid)

    def test_taxon_must_match_taxid(self, make_builder):
        with pytest.raises(GbkCreatorError, match="does not match"):
            make_builder(taxid=562)

    def test_email_required_for_lookup(self, monkeypatch):
        monkeypatch.delenv("NCBI_EMAIL", raising=False)
        with pytest.raises(TaxonomyError, match="NCBI_EMAIL"):
            GenBankBuilder(ID, 1308, strain="FAM13496")

    def test_invalid_transl_table(self, make_builder):
        with pytest.raises(GbkCreatorError, match="transl_table"):
            make_builder(transl_table=7)


class TestContigNames:
    def test_first_contig_must_be_scf1(self, make_builder):
        b = make_builder()
        with pytest.raises(NamingError, match="Expected 'FAM13496-i1-1.1_scf1'"):
            b.add_contig(f"{ID}_scf2", "ACGT")

    def test_contigs_must_be_numbered_in_add_order(self, make_builder):
        b = make_builder()
        b.add_contig(f"{ID}_scf1", "ACGT")
        with pytest.raises(NamingError):
            b.add_contig(f"{ID}_scf3", "ACGT")  # gap
        b.add_contig(f"{ID}_scf2", "ACGT")

    def test_duplicate_contig(self, make_builder):
        b = make_builder()
        b.add_contig(f"{ID}_scf1", "ACGT")
        with pytest.raises(NamingError, match="duplicate"):
            b.add_contig(f"{ID}_scf1", "ACGT")

    @pytest.mark.parametrize("name", [f"{ID}_scf01", f"{ID}_scf001", f"{ID}_scf0"])
    def test_no_leading_zeros(self, make_builder, name):
        with pytest.raises(NamingError, match="leading zeros"):
            make_builder().add_contig(name, "ACGT")

    @pytest.mark.parametrize("name", ["contig1", f"{ID}_contig1", f"{ID}scf1", "X_scf1", f"{ID}_scf1 "])
    def test_wrong_pattern(self, make_builder, name):
        with pytest.raises(NamingError):
            make_builder().add_contig(name, "ACGT")

    def test_chromosome_requires_circular(self, make_builder):
        b = make_builder()
        with pytest.raises(GbkCreatorError, match="must be circular"):
            b.add_contig(f"{ID}_scf1", "ACGT", topology="linear", kind="chromosome")
        b.add_contig(f"{ID}_scf1", "ACGT", topology="circular", kind="chromosome")

    def test_invalid_kind_and_topology(self, make_builder):
        b = make_builder()
        with pytest.raises(GbkCreatorError, match="kind"):
            b.add_contig(f"{ID}_scf1", "ACGT", kind="contig")
        with pytest.raises(GbkCreatorError, match="topology"):
            b.add_contig(f"{ID}_scf1", "ACGT", topology="ring")

    @pytest.mark.parametrize("seq", ["", "ACGT ACGT", "ACGU", "AC-GT", "ACGT\n"])
    def test_invalid_sequence(self, make_builder, seq):
        with pytest.raises(GbkCreatorError):
            make_builder().add_contig(f"{ID}_scf1", seq)

    def test_sequence_is_upper_cased_and_iupac_allowed(self, make_builder):
        b = make_builder()
        b.add_contig(f"{ID}_scf1", "acgtnrykm")
        assert b.contigs[0].sequence == "ACGTNRYKM"


class TestLocusTags:
    @pytest.mark.parametrize("lt", [f"{ID}_00001", f"{ID}_0000001", f"{ID}_00000a", f"{ID}000001",
                                    "OTHER_000001", f"{ID}_000001 ", f"{ID}_-00001"])
    def test_format(self, builder, lt):
        with pytest.raises(NamingError, match="6 digits"):
            builder.add_gene(lt, f"{ID}_scf1", loc(0, 30), "CDS", product="p", translation="M")

    def test_unique(self, builder):
        cds(builder, 1, 0, 30)
        with pytest.raises(NamingError, match="duplicate locus tag"):
            cds(builder, 1, 100, 130)

    def test_unknown_contig(self, builder):
        with pytest.raises(FeatureError, match="Unknown contig"):
            cds(builder, 1, 0, 30, contig=f"{ID}_scf9")

    def test_must_increase_along_position(self, builder):
        cds(builder, 2, 0, 30)
        cds(builder, 1, 100, 130)  # accepted when added; the order check is global
        with pytest.raises(NamingError, match="increase strictly"):
            builder.validate()

    def test_must_increase_along_contig_order(self, builder):
        cds(builder, 2, 0, 30, contig=f"{ID}_scf1")
        cds(builder, 1, 0, 30, contig=f"{ID}_scf2")
        with pytest.raises(NamingError, match="increase strictly"):
            builder.validate()

    def test_gaps_are_fine_and_add_order_is_free(self, builder, tmp_path):
        cds(builder, 5, 100, 130)
        cds(builder, 1, 0, 30)
        cds(builder, 900000, 0, 30, contig=f"{ID}_scf2")
        builder.validate()

    def test_wraparound_gene_sorts_first(self, make_builder):
        b = make_builder()
        b.add_contig(f"{ID}_scf1", random_dna(1000), "circular", "chromosome")
        cds(b, 2, 100, 400)
        wrap = CompoundLocation([loc(900, 1000, -1), loc(0, 50, -1)][::-1])
        b.add_gene(tag(1), f"{ID}_scf1", wrap, "CDS", product="p", translation="M")
        b.validate()


# ----------------------------------------------------------------------------- locations
class TestLocations:
    @pytest.mark.parametrize("strand", [None, 0])
    def test_strand_must_be_plus_or_minus(self, builder, strand):
        with pytest.raises(FeatureError, match="strand must be"):
            cds(builder, 1, 0, 30, strand=strand)

    def test_strand_of_every_join_part(self, builder):
        j = CompoundLocation([loc(0, 30, 1), SimpleLocation(40, 60, strand=None)])
        with pytest.raises(FeatureError, match="strand must be"):
            builder.add_gene(tag(1), f"{ID}_scf1", j, "CDS", product="p", translation="M")

    def test_join_parts_same_strand(self, builder):
        j = CompoundLocation([loc(0, 30, 1), loc(40, 60, -1)])
        with pytest.raises(FeatureError, match="same strand"):
            builder.add_gene(tag(1), f"{ID}_scf1", j, "CDS", product="p", translation="M")

    @pytest.mark.parametrize("start,end", [(0, 2001), (1990, 2100)])
    def test_must_lie_within_contig(self, builder, start, end):
        with pytest.raises(FeatureError, match="outside contig"):
            cds(builder, 1, start, end)

    def test_end_of_contig_ok(self, builder):
        cds(builder, 1, 1970, 2000)

    def test_must_be_biopython_location(self, builder):
        with pytest.raises(FeatureError, match="Biopython"):
            builder.add_gene(tag(1), f"{ID}_scf1", "1..30", "CDS", product="p", translation="M")

    def test_partial_ends_allowed(self, builder):
        part = SimpleLocation(BeforePosition(0), AfterPosition(30), strand=1)
        builder.add_gene(tag(1), f"{ID}_scf1", part, "CDS", product="p", translation="M", codon_start=2)

    def test_fuzzy_positions_rejected(self, builder):
        fuzzy = SimpleLocation(WithinPosition(5, left=5, right=8), 30, strand=1)
        with pytest.raises(FeatureError, match="unsupported start"):
            builder.add_gene(tag(1), f"{ID}_scf1", fuzzy, "CDS", product="p", translation="M")

    def test_only_join_operator(self, builder):
        order = CompoundLocation([loc(0, 30), loc(40, 60)], operator="order")
        with pytest.raises(FeatureError, match="only join"):
            builder.add_gene(tag(1), f"{ID}_scf1", order, "CDS", product="p", translation="M")
        with pytest.raises(FeatureError, match="only join"):
            builder.add_feature(f"{ID}_scf1", "misc_feature", order, note="x")

    def test_join_allowed(self, builder):
        j = CompoundLocation([loc(0, 30), loc(29, 60)])
        builder.add_gene(tag(1), f"{ID}_scf1", j, "CDS", product="p", translation="M",
                         ribosomal_slippage=True)


# ----------------------------------------------------------------------------- features
class TestFeatures:
    @pytest.mark.parametrize("ftype", ["gene", "misc_RNA", "mRNA", "exon", "cds", "source"])
    def test_gene_feature_types(self, builder, ftype):
        with pytest.raises(FeatureError, match="feature_type"):
            builder.add_gene(tag(1), f"{ID}_scf1", loc(0, 30), ftype, product="p")

    @pytest.mark.parametrize("ftype", ["CDS", "tRNA", "gene", "source", "mobile_element"])
    def test_untagged_feature_types(self, builder, ftype):
        with pytest.raises(FeatureError, match="feature_type"):
            builder.add_feature(f"{ID}_scf1", ftype, loc(0, 30))

    def test_untagged_types_ok(self, builder):
        builder.add_feature(f"{ID}_scf1", "repeat_region", loc(0, 30), rpt_family="CRISPR", rpt_type="direct",
                            rpt_unit_range="1..5", rpt_unit_seq="acgta", inference="x")
        builder.add_feature(f"{ID}_scf1", "regulatory", loc(40, 60), regulatory_class="riboswitch",
                            bound_moiety="adenosylcobalamin", db_xref="RFAM:RF00174", note="n")
        builder.add_feature(f"{ID}_scf1", "misc_feature", loc(70, 90), note="possible 16S")

    def test_untagged_has_no_locus_tag(self, builder):
        with pytest.raises(FeatureError, match="locus_tag"):
            builder.add_feature(f"{ID}_scf1", "repeat_region", loc(0, 30), locus_tag=tag(1))

    def test_regulatory_requires_class(self, builder):
        with pytest.raises(FeatureError, match="regulatory_class"):
            builder.add_feature(f"{ID}_scf1", "regulatory", loc(0, 30), note="x")

    def test_all_gene_types(self, builder):
        cds(builder, 1, 0, 30)
        cds(builder, 2, 40, 70, pseudo=True, translation=None)
        builder.add_gene(tag(3), f"{ID}_scf1", loc(100, 180), "tRNA", product="tRNA-Gly", anticodon="(pos:100..102,aa:Gly,seq:tcc)")
        builder.add_gene(tag(4), f"{ID}_scf1", loc(200, 300), "rRNA", product="16S ribosomal RNA")
        builder.add_gene(tag(5), f"{ID}_scf1", loc(400, 500), "ncRNA", product="RNase P", gene="rnpB", ncRNA_class="RNase_P_RNA", db_xref="RFAM:RF00010")
        builder.add_gene(tag(6), f"{ID}_scf1", loc(600, 700), "tmRNA", product="tmRNA", gene="ssrA", db_xref="RFAM:RF00023")
        builder.validate()

    def test_product_required(self, builder):
        for product in (None, "", "  ", "a\nb"):
            with pytest.raises(FeatureError, match="product"):
                builder.add_gene(tag(1), f"{ID}_scf1", loc(0, 30), "tRNA", product=product)

    def test_pseudo_cds_has_no_translation(self, builder):
        with pytest.raises(FeatureError, match="pseudo gene must not have a translation"):
            cds(builder, 1, 0, 30, pseudo=True)

    def test_non_pseudo_cds_requires_translation(self, builder):
        with pytest.raises(FeatureError, match="requires a translation"):
            cds(builder, 1, 0, 30, translation=None)

    @pytest.mark.parametrize("bad", ["False", "True", "", "yes", 0, 1, None])
    def test_pseudo_must_be_a_real_bool(self, builder, bad):
        with pytest.raises(FeatureError, match="pseudo must be True or False"):
            cds(builder, 1, 0, 30, pseudo=bad, translation=None)

    def test_pseudo_only_for_cds(self, builder):
        with pytest.raises(FeatureError, match="only CDS genes can be pseudo"):
            trna(builder, 1, 0, 30, pseudo=True)

    def test_translation_only_for_cds(self, builder):
        with pytest.raises(FeatureError, match="only CDS"):
            trna(builder, 1, 0, 30, translation="MKV")

    @pytest.mark.parametrize("bad", ["", "MK V", "MK1", "M\nK"])
    def test_translation_alphabet(self, builder, bad):
        with pytest.raises(FeatureError):
            cds(builder, 1, 0, 30, translation=bad)

    def test_translation_is_trusted_not_rechecked(self, builder):
        cds(builder, 1, 0, 30, translation="MKV" * 100)  # does not match the 30 nt: not our business

    def test_codon_start(self, builder):
        cds(builder, 1, 0, 30, codon_start=3)
        with pytest.raises(FeatureError, match="codon_start must be"):
            cds(builder, 2, 40, 70, codon_start=4)
        with pytest.raises(FeatureError, match="only valid for CDS"):
            trna(builder, 3, 100, 130, codon_start=1)

    def test_transl_table(self, make_builder):
        b = make_builder(transl_table=4)
        b.add_contig(f"{ID}_scf1", random_dna(100))
        cds(b, 1, 0, 30)
        assert b._genes[tag(1)].child_qualifiers["transl_table"] == [4]
        cds(b, 2, 40, 70, transl_table=11)
        with pytest.raises(FeatureError, match="transl_table"):
            cds(b, 3, 80, 99, transl_table=99)


class TestQualifiers:
    def test_ec_number_not_in_gbk(self, builder):
        with pytest.raises(FeatureError, match="EC numbers must not be written to the gbk"):
            cds(builder, 1, 0, 30, EC_number="1.1.1.1")

    @pytest.mark.parametrize("key", ["go_function", "go_process", "go_component"])
    def test_go_not_in_gbk(self, builder, key):
        with pytest.raises(FeatureError, match="GO terms must not be written"):
            cds(builder, 1, 0, 30, **{key: "GO:0000001"})

    def test_pseudo_via_argument_only(self, builder):
        with pytest.raises(FeatureError, match="pseudo=True argument of add_gene"):
            builder.add_feature(f"{ID}_scf1", "misc_feature", loc(0, 30), note="x", pseudo=True)

    @pytest.mark.parametrize("key", ["experiment", "function", "EC_number", "colour"])
    def test_only_oracle_qualifiers(self, builder, key):
        with pytest.raises(FeatureError, match="not allowed"):
            cds(builder, 1, 0, 30, **{key: "x"})

    def test_qualifier_scoping_per_feature_type(self, builder):
        with pytest.raises(FeatureError, match="anticodon"):
            cds(builder, 1, 0, 30, anticodon="x")
        with pytest.raises(FeatureError, match="protein_id"):
            trna(builder, 2, 100, 130, protein_id="x")
        with pytest.raises(FeatureError, match="ncRNA_class"):
            trna(builder, 3, 200, 230, ncRNA_class="x")

    @pytest.mark.parametrize("bad", ["", "GO:0000001", "RFAM:RF1", "taxon:1308", 5])
    def test_db_xref_is_rfam_only(self, builder, bad):
        with pytest.raises(FeatureError):
            builder.add_gene(tag(1), f"{ID}_scf1", loc(0, 30), "rRNA", product="p", db_xref=bad)

    def test_newlines_rejected(self, builder):
        with pytest.raises(FeatureError, match="control characters"):
            cds(builder, 1, 0, 30, note="a\nb")

    def test_multi_valued_and_flag(self, builder):
        cds(builder, 1, 0, 30, inference=["a", "b"], ribosomal_slippage=True)
        q = builder._genes[tag(1)].child_qualifiers
        assert q["inference"] == ["a", "b"] and q["ribosomal_slippage"] == [None]
        with pytest.raises(FeatureError, match="flag"):
            cds(builder, 2, 40, 70, ribosomal_slippage="yes")

    @pytest.mark.parametrize("gene", ["ab", "x" * 12, "a*b", "a.b", "a,b", "", "abc\n"])
    def test_gene_name_validated_against_gc_regex(self, builder, gene):
        with pytest.raises(FeatureError, match="GC regex"):
            cds(builder, 1, 0, 30, gene=gene)

    @pytest.mark.parametrize("gene", ["cas9", "abc", "x" * 11, "rpl_A-1", "a b'c", "dnaK/J"])
    def test_valid_gene_names(self, builder, gene):
        cds(builder, 1, 0, 30, gene=gene)

    def test_gene_validated_even_if_user_json_lacks_gc(self, make_builder, tmp_path):
        p = tmp_path / "cats.json"
        p.write_text('{"KG": "^K[0-9]{5}$"}')
        b = make_builder(annotation_regex=p)
        b.add_contig(f"{ID}_scf1", random_dna(100))
        with pytest.raises(FeatureError, match="GC regex"):
            cds(b, 1, 0, 30, gene="ab")

    def test_gene_parent_is_auto_created(self, builder):
        cds(builder, 1, 0, 30, gene="cas9", gene_synonym="x")
        g = builder._genes[tag(1)]
        assert g.gene_qualifiers == {"gene": ["cas9"], "locus_tag": [tag(1)], "gene_synonym": ["x"]}
        assert g.gene_location == g.location


# ----------------------------------------------------------------------------- annotations
class TestAnnotations:
    @pytest.fixture
    def b(self, builder):
        cds(builder, 1, 0, 30)
        cds(builder, 2, 40, 70, pseudo=True, translation=None)
        trna(builder, 3, 100, 130)
        return builder

    def test_unknown_locus_tag(self, b):
        with pytest.raises(AnnotationError, match="unknown locus tag"):
            b.add_annotation(tag(99), "KG", "K00001")

    def test_regex_mismatch(self, b):
        for anno_type, bad in (("KG", "K0001"), ("KG", "k00001"), ("EC", "1.1.1.1"), ("GO", "GO:123"),
                               ("KR", "K00001"), ("KG", "K00001\n")):
            with pytest.raises(AnnotationError, match="does not match|must not contain"):
                b.add_annotation(tag(1), anno_type, bad)

    def test_unknown_category_tells_user_to_add_to_json(self, b):
        with pytest.raises(AnnotationError, match="add it to your annotation JSON"):
            b.add_annotation(tag(1), "XX", "whatever")

    @pytest.mark.parametrize("anno_type,value", [("GC", "GC:cas9"), ("GP", "GP:something")])
    def test_gc_gp_forbidden(self, b, anno_type, value):
        with pytest.raises(AnnotationError, match="GC and GP sidecars are forbidden"):
            b.add_annotation(tag(1), anno_type, value)

    def test_cds_and_pseudo_cds_allowed(self, b):
        b.add_annotation(tag(1), "EC", "EC:1.1.1.1")
        b.add_annotation(tag(2), "EC", "EC:1.1.1.2")

    def test_rna_rejected(self, b):
        with pytest.raises(AnnotationError, match="only allowed on CDS genes"):
            b.add_annotation(tag(3), "KG", "K00001")

    def test_duplicates_ignored(self, b, tmp_path):
        for _ in range(3):
            b.add_annotation(tag(1), "KG", "K00001")
        b.write(tmp_path)
        assert (tmp_path / f"{ID}.KG").read_text() == f"{tag(1)}\tK00001\n"

    @pytest.mark.parametrize("bad", ["", "a,b", "a\tb", " x", "x ", 5, None])
    def test_unsafe_ids(self, builder, bad):
        cds(builder, 1, 0, 30)
        with pytest.raises(AnnotationError):
            builder.add_annotation(tag(1), "EP", bad)

    def test_custom_category_json_plain_and_arx_format(self, make_builder, tmp_path):
        plain = tmp_path / "plain.json"
        plain.write_text('{"XX": "^XX[0-9]+$"}')
        arx = tmp_path / "arx.json"
        arx.write_text('{"XX": {"name": "x", "color": "rgb(0,0,0)", "regex": "^XX[0-9]+$", "hyperlinks": []}}')
        for path in (plain, arx):
            b = make_builder(annotation_regex=path)
            b.add_contig(f"{ID}_scf1", random_dna(100))
            cds(b, 1, 0, 30)
            b.add_annotation(tag(1), "XX", "XX12")
            with pytest.raises(AnnotationError, match="does not match"):
                b.add_annotation(tag(1), "XX", "YY12")
            with pytest.raises(AnnotationError, match="unknown annotation category"):
                b.add_annotation(tag(1), "KG", "K00001")  # not in the user's JSON


# ----------------------------------------------------------------------------- global validation
class TestWrite:
    def test_needs_contig_and_gene(self, make_builder, tmp_path):
        b = make_builder()
        with pytest.raises(GbkCreatorError, match="No contigs"):
            b.write(tmp_path)
        b.add_contig(f"{ID}_scf1", random_dna(100))
        with pytest.raises(GbkCreatorError, match="No genes"):
            b.write(tmp_path)

    def test_nothing_is_written_if_invalid(self, builder, tmp_path):
        cds(builder, 2, 0, 30)
        cds(builder, 1, 100, 130)
        out = tmp_path / "out"
        with pytest.raises(NamingError):
            builder.write(out)
        assert not out.exists()

    def test_written_files_and_stale_sidecars(self, builder, tmp_path):
        cds(builder, 1, 0, 30)
        (tmp_path / f"{ID}.KG").write_text("stale\n")
        (tmp_path / f"{ID}.GC").write_text("stale\n")
        (tmp_path / "unrelated.KG").write_text("keep\n")
        paths = builder.write(tmp_path)
        assert [p.name for p in paths] == [f"{ID}.gbk", "organism.json", "genome.json"]
        assert sorted(p.name for p in tmp_path.iterdir()) == [f"{ID}.gbk", "genome.json", "organism.json", "unrelated.KG"]

    def test_creates_out_dir(self, builder, tmp_path):
        cds(builder, 1, 0, 30)
        builder.write(tmp_path / "a" / "b")
        assert (tmp_path / "a" / "b" / f"{ID}.gbk").is_file()

    @pytest.mark.parametrize("strain,identifier", [("my-plasmid-strain", "my-plasmid-strain-1")])
    def test_definition_must_not_change_arx_classification(self, make_builder, strain, identifier, tmp_path):
        b = make_builder(identifier, strain=strain)
        b.add_contig(f"{identifier}_scf1", random_dna(100), "linear", "scaffold")
        b.add_gene(f"{identifier}_000001", f"{identifier}_scf1", loc(0, 30), "CDS", product="p", translation="M")
        with pytest.raises(GbkCreatorError, match="classify"):
            b.write(tmp_path)


# ----------------------------------------------------------------------------- old_locus_tag
class TestOldLocusTag:
    def test_written_to_gene_and_child_as_last_qualifier(self, builder, tmp_path):
        from Bio import SeqIO
        cds(builder, 1, 0, 30, gene="cas9", old_locus_tag="OLD_0001")
        cds(builder, 2, 40, 70, pseudo=True, translation=None, old_locus_tag="OLD 0002/x")
        trna(builder, 3, 100, 130, old_locus_tag="OLD_0003")
        cds(builder, 4, 200, 230)
        builder.write(tmp_path)
        feats = next(SeqIO.parse(tmp_path / f"{ID}.gbk", "genbank")).features
        by = {(f.type, f.qualifiers.get("locus_tag", [None])[0]): f for f in feats}
        for n, old, ctype in ((1, "OLD_0001", "CDS"), (2, "OLD 0002/x", "CDS"), (3, "OLD_0003", "tRNA")):
            for ftype in ("gene", ctype):
                q = by[ftype, tag(n)].qualifiers
                assert q["old_locus_tag"] == [old]
                assert list(q)[-1] == "old_locus_tag"  # after /translation too, like arx-tools' rename mode
        assert "old_locus_tag" not in by["gene", tag(4)].qualifiers and "old_locus_tag" not in by["CDS", tag(4)].qualifiers
        assert list(by["CDS", tag(1)].qualifiers)[-2:] == ["translation", "old_locus_tag"]

    @pytest.mark.parametrize("bad", ["", " ", " x", "x ", 'a"b', "a\nb", "a\tb", "x" * 101, 5, ["x"], b"x"])
    def test_validation(self, builder, bad):
        with pytest.raises(NamingError, match="old_locus_tag"):
            cds(builder, 1, 0, 30, old_locus_tag=bad)
        assert tag(1) not in builder._genes  # nothing half-added

    def test_max_length_100(self, builder):
        cds(builder, 1, 0, 30, old_locus_tag="o" * 100)

    def test_must_differ_from_new_locus_tag(self, builder):
        with pytest.raises(NamingError, match="must differ"):
            cds(builder, 1, 0, 30, old_locus_tag=tag(1))

    def test_unique_across_the_genome(self, builder):
        cds(builder, 1, 0, 30, old_locus_tag="OLD_1")
        with pytest.raises(NamingError, match=rf"already used by {tag(1)}"):
            cds(builder, 2, 40, 70, old_locus_tag="OLD_1", contig=f"{ID}_scf2")
        assert tag(2) not in builder._genes
        cds(builder, 2, 40, 70, old_locus_tag="OLD_2")

    def test_failed_add_does_not_reserve_the_old_tag(self, builder):
        with pytest.raises(FeatureError):
            cds(builder, 1, 0, 30, old_locus_tag="OLD_1", translation=None)  # non-pseudo CDS without translation
        cds(builder, 1, 0, 30, old_locus_tag="OLD_1")

    def test_extra_qualifier_hints_to_the_argument(self, builder):
        with pytest.raises(FeatureError, match="old_locus_tag argument of add_gene"):
            builder.add_feature(f"{ID}_scf1", "misc_feature", loc(0, 30), note="x", old_locus_tag="OLD")
        with pytest.raises(FeatureError, match="old_locus_tag argument of add_gene"):
            builder._normalize_qualifiers("CDS", {"old_locus_tag": ["x"]}, set(), "x")
