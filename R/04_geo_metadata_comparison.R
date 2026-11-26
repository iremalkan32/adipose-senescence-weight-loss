# GSE59034 ve GSE199063 GEO kayıtlarının ayrıntılı karşılaştırması (yalnızca metadata; birkaç dakika sürer)
options(timeout = 600, width = 250)
indir <- function(acc) {
  kok <- paste0(substr(acc, 1, nchar(acc) - 3), "nnn")
  hedef <- paste0(acc, "_series_matrix.txt.gz")
  for (u in c(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", kok, "/", acc, "/matrix/", hedef),
              paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", kok, "/", acc, "/matrix/", hedef)))
    if (!file.exists(hedef) || file.size(hedef) == 0) try(download.file(u, hedef, mode = "wb", quiet = TRUE), silent = TRUE)
  readLines(gzfile(hedef), n = 2000, warn = FALSE)
}
oku <- function(acc) {
  s <- indir(acc); s <- s[seq_len(which(grepl("^!series_matrix_table_begin", s))[1] - 1)]
  seri <- grep("^!Series_", s, value = TRUE)
  seri <- seri[!grepl("^!Series_(sample_id|supplementary_file|platform_taxid|sample_taxid)", seri)]
  ornek <- grep("^!Sample_", s, value = TRUE)
  alan <- sub("\t.*", "", ornek)
  deger <- lapply(strsplit(ornek, "\t"), function(v) gsub('^"|"$', "", v[-1]))
  # Aynı adlı alanlar (characteristics_ch1 gibi) numaralandırılır
  alan <- make.unique(sub("^!Sample_", "", alan), sep = ".")
  tab <- as.data.frame(setNames(deger, alan), check.names = FALSE, stringsAsFactors = FALSE)
  list(seri = seri, tab = tab)
}
A <- oku("GSE59034"); B <- oku("GSE199063")

for (x in list(list("GSE59034", A), list("GSE199063", B))) {
  cat("\n\n######################## ", x[[1]], " — SERİ BİLGİLERİ ########################\n")
  cat(x[[2]]$seri, sep = "\n")
}
gostergeler <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "source_name_ch1",
                 grep("^characteristics_ch1", colnames(A$tab), value = TRUE), "description", "contact_name",
                 "contact_institute", "contact_department", "contact_laboratory", "platform_id", "scan_protocol", "extract_protocol_ch1",
                 "data_processing")
for (x in list(list("GSE59034", A$tab), list("GSE199063", B$tab))) {
  t <- x[[2]]
  cat("\n\n######################## ", x[[1]], " — ÖRNEK ALANLARI (sütun adı: farklı değer sayısı) ########################\n")
  for (c in colnames(t)) {
    u <- unique(t[[c]])
    cat(sprintf("%-28s [%3d farklı] %s\n", c, length(u), substr(paste(head(u, 6), collapse = " | "), 1, 220)))
  }
}
cat("\n\n######################## GSE59034 — TÜM ÖRNEKLER ########################\n")
sut <- intersect(c("title", "geo_accession", "submission_date", "source_name_ch1",
                   grep("^characteristics_ch1", colnames(A$tab), value = TRUE), "description"), colnames(A$tab))
print(A$tab[, sut], row.names = FALSE)
cat("\n\n######################## GSE199063 — TÜM ÖRNEKLER ########################\n")
sut <- intersect(c("title", "geo_accession", "submission_date", "source_name_ch1",
                   grep("^characteristics_ch1", colnames(B$tab), value = TRUE), "description"), colnames(B$tab))
print(B$tab[, sut], row.names = FALSE)
write.csv(A$tab, "GSE59034_ornek_bilgileri.csv", row.names = FALSE)
write.csv(B$tab, "GSE199063_ornek_bilgileri.csv", row.names = FALSE)
cat("\nBitti.\n")
