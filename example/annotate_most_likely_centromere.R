# This script annotates single breakends to centromeres with the most likely centromeric breakends
# Requires a RepeatMasker annotated GRIDSS VCF
library(StructuralVariantAnnotation)
library(tidyverse)
library(argparser)
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="Annotated VCF")
argp = add_argument(argp, "--ref", default="", help="Reference genome to use. hg38 or a T2T complete genome is recommended")
# argv = parse_args(argp, argv=c("--ref", "../../ref/chm13.draft_v1.0.fasta", "--input", "DO52605T.purple.sv.vcf", "--output", "out.DO52605T.centromeric.vcf"))
argv = parse_args(argp)

vcf = readVcf(argv$input)
begr = breakendRanges(vcf)
centrogr = begr[info(vcf[begr$sourceId])$INSRMRC %in% c("Satellite/centr") | info(vcf[begr$sourceId])$INSRMRT %in% c("HSATII", "(CATTC)n", "(GAATG)n")]

fa = DNAStringSet(str_replace(centrogr$ALT, stringr::fixed("."), ""))
names(fa) = names(centrogr)
tmp_in  = paste0(argv$input, ".tmp.blatin.fasta")
tmp_out = paste0(argv$input, ".tmp.blatout.psl")
export(fa, tmp_in)
cmd = paste0("blat -fastMap -noHead ", argv$ref, " ", tmp_in, " ", tmp_out)
write(paste("Running", cmd), stderr())
system(cmd)
if (!file.exists(tmp_out)) {
	stop(paste(tmp_out, " does not exist. blat failed."))
}
alndf = read_tsv(
		tmp_out,
		col_names=c("match", "mismatch", "repmatch", "Ns", "Qcount", "Qbases", "Tcount", "Tbases", "strand", "Qname", "Qsize", "Qstart", "Qend", "Tname", "Tsize", "Tstart", "Tend", "blockcount", "blockSizes", "Qstarts", "tStarts"),
	col_types=c("iiiiiiiicciiiciiiiccc")) %>%
	mutate(
		QSizeMult=9-floor(Qsize/100),
		score=(1000-(QSizeMult*mismatch+Qcount+Tcount))*pmin(match/Qsize,1))
best_alndf = alndf %>%
	filter(!str_detect(Tname, stringr::fixed("chrUn")) & !str_detect(Tname, stringr::fixed("_random"))) %>%
	group_by(Qname, Tname) %>%
	arrange(desc(score)) %>%
	slice_head() %>%
	group_by(Qname) %>%
	arrange(desc(score)) %>%
	mutate(nextChrScore=lead(score)) %>%
	slice_head()
info(header(vcf)) = unique(as(rbind(
		as.data.frame(info(header(vcf))),
		data.frame(row.names="MLCENTRO", Number=1, Type="String",Description="Most likely centromere")),
	"DataFrame"))
best_alndf = best_alndf %>% 
	filter(score >= nextChrScore + 25)
info(vcf)$MLCENTRO = NA_character_
info(vcf[best_alndf$Qname])$MLCENTRO = best_alndf$Tname
writeVcf(vcf, argv$output)
