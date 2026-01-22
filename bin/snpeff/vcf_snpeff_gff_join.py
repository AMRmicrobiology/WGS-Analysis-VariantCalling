#!/usr/bin/env python3
import argparse
import re
import sys
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Une VCF SnpEff con GFF3 Prokka/Bakta + contexto intergénico mejorado (B, C y D)")
    parser.add_argument("--vcf", required=True, help="Archivo VCF anotado")
    parser.add_argument("--gff", required=True, help="Archivo GFF3 enriquecido")
    parser.add_argument("--out", required=True, help="Nombre del archivo TSV de salida")
    return parser.parse_args()


def load_gff(gff_path):
    """Carga el GFF en memoria con priorización de Name > product > prokka_product."""
    genes = []

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attrs_str = cols
            if feature not in ("CDS", "gene"):
                continue

            attrs = {}
            for a in attrs_str.split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k] = v

            # PRIORIDAD: Name > product > prokka_product
            best_product = attrs.get("Name", attrs.get("product", attrs.get("prokka_product", ".")))
            
            # PRIORIDAD: gene > prokka_gene
            best_gene = attrs.get("gene", attrs.get("prokka_gene", "."))

            entry = {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "ID": attrs.get("ID", ""),
                "locus_tag": attrs.get("locus_tag", ""),
                "gene": best_gene,
                "product": best_product,
                "prokka_gene": attrs.get("prokka_gene", "."),
                "prokka_product": attrs.get("prokka_product", ".")
            }

            genes.append(entry)

    return genes


def find_neighbors(genes, chrom, pos):
    """Devuelve (closest, left, right) + distancias."""
    pos = int(pos)
    same_chrom = [g for g in genes if g["chrom"] == chrom]

    left = None
    right = None

    for g in same_chrom:
        if g["end"] < pos:
            if left is None or g["end"] > left["end"]:
                left = g
        if g["start"] > pos:
            if right is None or g["start"] < right["start"]:
                right = g

    # Gen más cercano
    closest = None
    dist_closest = None

    if left:
        d_left = pos - left["end"]
    else:
        d_left = None

    if right:
        d_right = right["start"] - pos
    else:
        d_right = None

    if d_left is not None and (d_right is None or d_left <= d_right):
        closest = left
        dist_closest = d_left
    elif d_right is not None:
        closest = right
        dist_closest = d_right

    return closest, left, right, d_left, d_right


def find_overlapping_gene(genes, chrom, pos):
    """Busca si la variante cae DENTRO de un gen."""
    pos = int(pos)
    for g in genes:
        if g["chrom"] == chrom and g["start"] <= pos <= g["end"]:
            return g
    return None


def main():
    args = parse_args()

    # ---------------------------
    # 1. CARGAR GFF
    # ---------------------------
    if not os.path.exists(args.gff):
        sys.exit(f"ERROR: GFF no encontrado: {args.gff}")

    genes = load_gff(args.gff)
    print(f"[INFO] Cargados {len(genes)} genes del GFF", file=sys.stderr)

    # ---------------------------
    # 2. PREPARAR ARCHIVOS DE SALIDA
    # ---------------------------
    if not os.path.exists(args.vcf):
        sys.exit(f"ERROR: VCF no encontrado: {args.vcf}")

    # Generar nombre del reporte humano
    base_name = args.out.replace(".tsv", "").replace(".txt", "")
    report_file = f"{base_name}.report.txt"

    ann_regex = re.compile(r"ANN=([^;]+)")

    header = [
        "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER",
        "Gene_Name", "Gene_ID", "IMPACT", "Effect", "HGVSc", "HGVSp",
        "gene_name", "gene_product", "locus_tag",
        "closest_gene", "closest_product", "closest_locus", "distance_closest",
        "left_gene", "left_product", "left_locus", "left_distance",
        "right_gene", "right_product", "right_locus", "right_distance",
        "flanking_genes", "flanking_products"
    ]

    # Abrir ambos archivos
    out_tsv = open(args.out, "w")
    out_report = open(report_file, "w")

    # Header TSV
    out_tsv.write("\t".join(header) + "\n")

    # Header Reporte Humano
    out_report.write("=" * 80 + "\n")
    out_report.write("REPORTE DE ANOTACIÓN DE VARIANTES\n")
    out_report.write(f"VCF: {os.path.basename(args.vcf)}\n")
    out_report.write(f"GFF: {os.path.basename(args.gff)}\n")
    out_report.write("=" * 80 + "\n\n")

    # ---------------------------
    # 3. PROCESAR VCF
    # ---------------------------
    with open(args.vcf) as vcf:
        for line in vcf:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt, qual, fil, info = cols[:8]

            # ============================
            # PARSEAR ANN DE SNPEFF
            # ============================
            m = ann_regex.search(info)
            if m:
                ann_raw = m.group(1)
                first = ann_raw.split(",")[0]
                parts = first.split("|")

                gene_name = parts[3] if len(parts) > 3 else "."
                gene_id = parts[4] if len(parts) > 4 else "."
                effect = parts[1] if len(parts) > 1 else "."
                impact = parts[2] if len(parts) > 2 else "."
                hgvsc = parts[9] if len(parts) > 9 else "."
                hgvsp = parts[10] if len(parts) > 10 else "."
            else:
                gene_name = gene_id = effect = impact = hgvsc = hgvsp = "."

            # ============================
            # BUSCAR GEN AFECTADO (MATCH DIRECTO)
            # ============================
            # Primero buscar por ID exacto
            match = None
            for g in genes:
                if gene_id == g["ID"] or gene_name == g["locus_tag"] or gene_name == g["gene"]:
                    match = g
                    break
            
            # Si no hay match por ID, buscar si la variante cae dentro de un gen
            if not match:
                match = find_overlapping_gene(genes, chrom, pos)

            if match:
                p_gene = match["gene"]
                p_prod = match["product"]
                p_locus = match["locus_tag"]
            else:
                p_gene = p_prod = p_locus = "."

            # ============================
            # OPCIONES B + C + D (VECINOS)
            # ============================
            closest, left, right, d_left, d_right = find_neighbors(genes, chrom, pos)

            # closest
            if closest:
                closest_gene = closest["gene"]
                closest_product = closest["product"]
                closest_locus = closest["locus_tag"]
                dist_closest = d_left if left == closest else d_right
            else:
                closest_gene = closest_product = closest_locus = dist_closest = "."

            # left
            if left:
                left_gene = left["gene"]
                left_prod = left["product"]
                left_loc = left["locus_tag"]
                left_dist = d_left
            else:
                left_gene = left_prod = left_loc = left_dist = "."

            # right
            if right:
                right_gene = right["gene"]
                right_prod = right["product"]
                right_loc = right["locus_tag"]
                right_dist = d_right
            else:
                right_gene = right_prod = right_loc = right_dist = "."

            # flanking lists
            flank_genes = f"{left_gene}|{right_gene}" if left_gene != "." or right_gene != "." else "."
            flank_products = f"{left_prod}|{right_prod}" if left_prod != "." or right_prod != "." else "."

            # ============================
            # ESCRIBIR FILA TSV
            # ============================
            row = [
                chrom, pos, ref, alt, qual, fil,
                gene_name, gene_id, impact, effect, hgvsc, hgvsp,
                p_gene, p_prod, p_locus,
                closest_gene, closest_product, closest_locus, str(dist_closest),
                left_gene, left_prod, left_loc, str(left_dist),
                right_gene, right_prod, right_loc, str(right_dist),
                flank_genes, flank_products
            ]

            out_tsv.write("\t".join(str(x) for x in row) + "\n")

            # ============================
            # ESCRIBIR BLOQUE HUMANO
            # ============================
            out_report.write(f"▶ VARIANTE: {chrom}:{pos} {ref}→{alt}\n")
            out_report.write(f"  Calidad: {qual} | Filtro: {fil}\n")
            out_report.write(f"  Impacto: {impact} | Efecto: {effect}\n")
            out_report.write(f"  HGVS.c: {hgvsc}\n")
            out_report.write(f"  HGVS.p: {hgvsp}\n")
            out_report.write(f"  SnpEff Gene: {gene_name} ({gene_id})\n")
            
            # Si hay match directo, mostrar info del gen afectado
            if match:
                out_report.write(f"  ✓ GEN AFECTADO: {p_gene} | {p_prod} | {p_locus}\n")
            
            out_report.write("-" * 80 + "\n")

            out_report.write(f"closest_gene_name:         {closest_gene}\n")
            out_report.write(f"closest_gene_product:      {closest_product}\n")
            out_report.write(f"closest_locus_tag:         {closest_locus}\n")
            out_report.write(f"distance_to_closest:       {dist_closest} bp\n\n")

            out_report.write(f"left_gene_name:            {left_gene}\n")
            out_report.write(f"left_gene_product:         {left_prod}\n")
            out_report.write(f"left_locus_tag:            {left_loc}\n")
            out_report.write(f"left_distance:             {left_dist} bp\n\n")

            out_report.write(f"right_gene_name:           {right_gene}\n")
            out_report.write(f"right_gene_product:        {right_prod}\n")
            out_report.write(f"right_locus_tag:           {right_loc}\n")
            out_report.write(f"right_distance:            {right_dist} bp\n\n")

            out_report.write(f"genes_flanking:            {flank_genes}\n")
            out_report.write(f"products_flanking:         {flank_products}\n")
            out_report.write("=" * 80 + "\n\n")

    out_tsv.close()
    out_report.close()

    print(f"[INFO] TSV generado → {args.out}", file=sys.stderr)
    print(f"[INFO] Reporte humano generado → {report_file}", file=sys.stderr)


if __name__ == "__main__":
    main()