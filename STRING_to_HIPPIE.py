import pandas as pd


GENE_INFO_COLS = [
    "tax_id",
    "GeneID",
    "Symbol",
    "LocusTag",
    "Synonyms",
    "dbXrefs",
    "chromosome",
    "map_location",
    "description",
    "type_of_gene",
    "Symbol_from_nomenclature_authority",
    "Full_name_from_nomenclature_authority",
    "Nomenclature_status",
    "Other_designations",
    "Modification_date",
    "Feature_type",
]


def ordered_unique(seq):
    seen = set()
    out = []
    for x in seq:
        if x and x != "-" and x not in seen:
            seen.add(x)
            out.append(x)

    return out


def load_seed_ids(path):
    seed_ids = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                seed_ids.append(parts[0])

    return ordered_unique(seed_ids)


def load_string_aliases(path):
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["protein_id", "alias", "source"],
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        compression="infer"
    )

    df["protein_id"] = df["protein_id"].str.strip()
    df["alias"] = df["alias"].str.strip()
    df["source"] = df["source"].str.strip()

    return df


def load_gene_info(path):
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=GENE_INFO_COLS,
        dtype=str,
        keep_default_na=False,
        low_memory=False
    )

    df = df[df["tax_id"] == "9606"].copy()

    df["canonical_symbol"] = df["Symbol"]
    mask = (
        df["Symbol_from_nomenclature_authority"].ne("") &
        df["Symbol_from_nomenclature_authority"].ne("-")
    )

    df.loc[mask, "canonical_symbol"] = df.loc[
        mask, "Symbol_from_nomenclature_authority"
    ]
    
    return df


def build_gene_maps(gene_info_df):
    geneid_map = {}
    official_map = {}
    synonym_map = {}

    for row in gene_info_df.itertuples(index=False):
        gene_id = row.GeneID
        canonical_symbol = row.canonical_symbol

        geneid_map[gene_id] = (gene_id, canonical_symbol)

        symbol_candidates = ordered_unique([
            row.Symbol,
            row.Symbol_from_nomenclature_authority,
            canonical_symbol,
        ])

        for sym in symbol_candidates:
            official_map[sym.upper()] = (gene_id, canonical_symbol)

        if row.Synonyms and row.Synonyms != "-":
            for syn in row.Synonyms.split("|"):
                syn = syn.strip()
                if syn and syn != "-":
                    synonym_map.setdefault(syn.upper(), []).append(
                        (gene_id, canonical_symbol)
                    )

    for syn, vals in synonym_map.items():
        synonym_map[syn] = ordered_unique(vals)

    return geneid_map, official_map, synonym_map


def resolve_aliases(aliases, geneid_map, official_map, synonym_map):
    aliases = ordered_unique(aliases)

    for alias in aliases:
        if alias in geneid_map:
            return geneid_map[alias], "geneid_alias", alias

    for alias in aliases:
        hit = official_map.get(alias.upper())
        if hit is not None:
            return hit, "symbol_alias", alias

    for alias in aliases:
        hits = synonym_map.get(alias.upper(), [])
        if len(hits) == 1:
            return hits[0], "unique_synonym_alias", alias

    return (None, None), "unmapped", None


def load_hippie_nodes(hippie_mitab_path):
    df = pd.read_csv(
        hippie_mitab_path,
        sep="\t",
        header=0,
        usecols=[0, 1, 16, 17],
        dtype=str,
        keep_default_na=False,
        low_memory=False
    )

    df.columns = ["id_a", "id_b", "symbol_a", "symbol_b"]

    entrez_a = df["id_a"].str.extract(r"^entrez gene:(\d+)$", expand=False)
    entrez_b = df["id_b"].str.extract(r"^entrez gene:(\d+)$", expand=False)

    valid_entrez = set(
        pd.concat([entrez_a, entrez_b], ignore_index=True).dropna().tolist()
    )

    valid_symbols = set(
        pd.concat([df["symbol_a"], df["symbol_b"]], ignore_index=True)
        .replace("-", pd.NA)
        .dropna()
        .str.upper()
        .tolist()
    )

    return valid_entrez, valid_symbols


def main():
    # ../9606.protein.aliases.v12.0.txt.gz
    string_alias_file = input("Enter STRING alias file path/name: ").strip()

    # data/seed_nodes/genetic/parkinsons_disease.txt
    seed_file = input("Enter STRING seed file path/name: ").strip()

    # ../Homo_sapiens.gene_info
    gene_info_file = input("Enter path to Homo_sapiens.gene_info: ").strip()

    # ../HIPPIE_network.txt
    hippie_mitab_file = input("Enter path to HIPPIE current MITAB: ").strip()

    seed_ids = load_seed_ids(seed_file)
    alias_df = load_string_aliases(string_alias_file)
    gene_info_df = load_gene_info(gene_info_file)
    valid_hippie_entrez, valid_hippie_symbols = load_hippie_nodes(hippie_mitab_file)

    geneid_map, official_map, synonym_map = build_gene_maps(gene_info_df)

    alias_groups = (
        alias_df[alias_df["protein_id"].isin(seed_ids)]
        .groupby("protein_id")["alias"]
        .apply(list)
        .to_dict()
    )

    rows = []

    for pid in seed_ids:
        aliases = alias_groups.get(pid, [])
        (gene_id, symbol), match_type, matched_alias = resolve_aliases(
            aliases,
            geneid_map,
            official_map,
            synonym_map
        )

        hippie_entrez_label = ""
        present_in_hippie = False

        if gene_id is not None:
            hippie_entrez_label = f"entrez gene:{gene_id}"
            present_in_hippie = (
                gene_id in valid_hippie_entrez or
                symbol.upper() in valid_hippie_symbols
            )

        rows.append({
            "string_protein_id": pid,
            "matched_alias": matched_alias if matched_alias else "",
            "match_type": match_type,
            "hippie_gene_id": gene_id if gene_id else "",
            "hippie_entrez_label": hippie_entrez_label,
            "hippie_symbol": symbol if symbol else "",
            "present_in_hippie": present_in_hippie,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv("string_to_hippie_mapping.tsv", sep="\t", index=False)

    present_df = out_df[
        (out_df["present_in_hippie"]) &
        (out_df["hippie_gene_id"] != "")
    ].copy()

    present_df["hippie_entrez_label"].to_csv(
        "hippie_seed_entrez.txt",
        index=False,
        header=False
    )

    present_df["hippie_symbol"].to_csv(
        "hippie_seed_symbols.txt",
        index=False,
        header=False
    )

    missing_df = out_df[~out_df["present_in_hippie"]].copy()
    missing_df.to_csv("string_to_hippie_unmatched.tsv", sep="\t", index=False)

    print("Saved:")
    print("  string_to_hippie_mapping.tsv")
    print("  hippie_seed_entrez.txt")
    print("  hippie_seed_symbols.txt")
    print("  string_to_hippie_unmatched.tsv")


if __name__ == "__main__":
    main()