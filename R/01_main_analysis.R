###############################################################################
# PİLOT ANALİZ: Kilo kaybı sonrası yağ dokusu senesans imzası
# Veri setleri: GSE199063 (bariatrik cerrahi; öncesi, 2. yıl, 5. yıl)
#               GSE77962  (diyet; öncesi, kilo kaybı sonrası, idame sonrası)
#
# Soru: SenMayo senesans skoru kilo kaybından sonra düşüyor mu?
#       Düşüş, rastgele gen setlerinden beklenenden büyük mü?
#
# Kullanım: R veya RStudio'da bu dosyayı açın ve tamamını çalıştırın
#           (RStudio: Ctrl+Shift+Enter / Cmd+Shift+Enter).
# Tüm çıktılar "senesans_pilot_cikti" klasörüne kaydedilir.
###############################################################################

## ---------------------------------------------------------------------------
## 0) AYARLAR
## ---------------------------------------------------------------------------
CIKTI_KLASORU <- "senesans_pilot_cikti"
VERI_SETLERI  <- c("GSE199063", "GSE77962", "GSE59034", "GSE65540")
RNASEQ_SETLERI <- c("GSE65540")   # ifade verisi NCBI'nin hazırladığı sayım matrisinden okunur
RASTGELE_SET_SAYISI <- 1000   # özgüllük testi için rastgele gen seti sayısı
DUYARLILIK <- TRUE            # HPA eşik ve tek gen çıkarma duyarlılık analizleri
PROB_BIRLESTIRME <- "en_yuksek"  # aynı gene ait problar: "en_yuksek" (ortalaması en yüksek prob) veya "ortalama"
set.seed(2026)

# Otomatik algılama başarısız olursa buraya sütun adlarını elle yazın.
# Sütun adlarını "01_metadata_<GSE>.txt" dosyasında görebilirsiniz.
# Örnek: ELLE_AYAR$GSE199063 <- list(zaman = "time point:ch1", kisi = "subject:ch1",
#                                    baslangic = "before")
# Makrofaj işaret genleri (SenMayo ile ortak gen içermez)
MAKROFAJ_GENLERI <- c("CD68", "CD14", "CD163", "ITGAM", "ITGAX", "TREM2", "LYVE1", "MRC1", "MSR1",
                      "C1QA", "C1QB", "C1QC", "CSF1R", "FCGR3A", "LIPA", "PLA2G7", "ADGRE1", "MS4A7",
                      "TYROBP", "FCER1G", "LAPTM5", "AIF1", "CD86", "SIGLEC1")
T_HUCRE_GENLERI <- c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "CD8A", "CD8B", "CD247", "LCK", "ZAP70", "GZMK", "TRAT1", "CD96")
MAST_GENLERI    <- c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC", "HPGDS", "CMA1", "SIGLEC6", "CLC")
ADIPOSIT_GENLERI <- c("ADIPOQ", "PLIN1", "PLIN4", "LEP", "CIDEA", "CIDEC", "LIPE", "PNPLA2", "FABP4", "GPD1", "TUSC5", "PPARG")
MAKROFAJ_KOR_ESIGI <- 0.5   # bir SenMayo geni herhangi bir bağışıklık skoruyla bu düzeyin üstünde ilişkiliyse "immün" sayılır   # bu korelasyonun üstündeki SenMayo genleri "makrofajla ilişkili" sayılır
# Klasik senesans belirteçleri (LMNB1 senesansta azalır)
TEMEL_SENESANS_GENLERI <- c("CDKN1A", "CDKN2A", "CDKN2B", "TP53", "LMNB1", "GLB1", "SERPINE1", "MKI67")

ELLE_AYAR <- list(
  GSE199063 = list(zaman = "time_point:ch1", kisi = "subject_id:ch1",
                   baslangic = "Baseline", kontrol = "Control"),
  GSE77962  = list(zaman = "time point:ch1", kisi = "subject_id:ch1",
                   baslangic = "at study start", grup = "treatment:ch1")
)

## ---------------------------------------------------------------------------
## 1) PAKETLER
## ---------------------------------------------------------------------------
# Gen seti dosyaları (senmayo.txt, fridman_up.txt) klasörde varsa msigdbr gerekmez
DOSYADAN_GEN_SETI <- file.exists("senmayo.txt")
cran <- c("ggplot2", if (!DOSYADAN_GEN_SETI) "msigdbr")
bioc <- c("GEOquery", "limma", "Biobase")
eksik_cran <- cran[!sapply(cran, requireNamespace, quietly = TRUE)]
if (length(eksik_cran)) install.packages(eksik_cran)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
eksik_bioc <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(eksik_bioc)) BiocManager::install(eksik_bioc, update = FALSE, ask = FALSE)
# msigdbr 10+ sürümünde gen setleri ayrı "msigdbdf" paketinde
if (!DOSYADAN_GEN_SETI && packageVersion("msigdbr") >= "10.0.0" && !requireNamespace("msigdbdf", quietly = TRUE)) {
  install.packages("msigdbdf",
                   repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))
}

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(Biobase); library(ggplot2)
})
options(timeout = 3600)  # büyük dosyalar için indirme süresi
dir.create(CIKTI_KLASORU, showWarnings = FALSE)
yol <- function(...) file.path(CIKTI_KLASORU, paste0(...))

## ---------------------------------------------------------------------------
## 2) GEN SETLERİ (MSigDB)
## ---------------------------------------------------------------------------
gen_setlerini_al <- function() {
  oku <- function(f) if (file.exists(f)) unique(trimws(readLines(f))) else NULL
  if (DOSYADAN_GEN_SETI) {
    setler <- list(SenMayo = oku("senmayo.txt"), Fridman_SEN_UP = oku("fridman_up.txt"))
    return(setler[!sapply(setler, is.null)])
  }
  library(msigdbr)
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "C2"))
  bul <- function(desen) {
    ad <- unique(m$gs_name[grepl(desen, m$gs_name)])
    if (!length(ad)) return(NULL)
    unique(m$gene_symbol[m$gs_name == ad[1]])
  }
  setler <- list(SenMayo = bul("SEN_?MAYO"), Fridman_SEN_UP = bul("^FRIDMAN_SENESCENCE_UP$"))
  setler <- setler[!sapply(setler, is.null)]
  if (is.null(setler$SenMayo))
    stop("SenMayo gen seti bulunamadı. Klasöre senmayo.txt dosyası koyun.")
  setler
}

## ---------------------------------------------------------------------------
## 3) VERİ İNDİRME VE GEN SEMBOLÜNE ÇEVİRME
## ---------------------------------------------------------------------------
sembol_bul <- function(fd) {
  aday <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "gene_symbol",
            "GeneSymbol")
  s <- intersect(aday, colnames(fd))
  if (length(s)) return(trimws(sub(" ///.*", "", as.character(fd[[s[1]]]))))
  if ("gene_assignment" %in% colnames(fd)) {         # Affymetrix Gene/Clariom dizileri
    return(trimws(vapply(strsplit(as.character(fd$gene_assignment), " // "),
                         function(x) if (length(x) >= 2) x[2] else NA_character_, "")))
  }
  NULL
}

# GEO web sayfası (www.ncbi.nlm.nih.gov) otomatik isteklere doğrulama ekranı gösterebildiği için
# dosyalar doğrudan NCBI dosya sunucusundan (ftp.ncbi.nlm.nih.gov) indirilir.
ftp_klasor <- function(id, tur) {             # GSE77962 -> .../series/GSE77nnn/GSE77962/
  n <- sub("^[A-Z]+", "", id)
  kok <- if (nchar(n) > 3) substr(n, 1, nchar(n) - 3) else ""
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/%s/%s%snnn/%s/", tur, substr(id, 1, 3), kok, id)
}
INDIRME_KLASORU <- "geo_indirilen"
dir.create(INDIRME_KLASORU, showWarnings = FALSE)
indir <- function(url) {
  hedef <- file.path(INDIRME_KLASORU, basename(url))
  if (file.exists(hedef) && file.size(hedef) > 0) return(hedef)
  for (u in c(url, sub("^https://", "ftp://", url))) {
    ok <- tryCatch({ download.file(u, hedef, mode = "wb", quiet = TRUE); file.size(hedef) > 0 },
                   error = function(e) { message("    indirilemedi: ", u, " (", conditionMessage(e), ")"); FALSE },
                   warning = function(w) { message("    indirilemedi: ", u, " (", conditionMessage(w), ")"); FALSE })
    if (isTRUE(ok)) return(hedef)
  }
  if (file.exists(hedef)) file.remove(hedef)
  NULL
}
dizin_listesi <- function(url, desen) {
  html <- tryCatch(paste(readLines(url, warn = FALSE), collapse = ""), error = function(e) "")
  unique(regmatches(html, gregexpr(desen, html))[[1]])
}
# Platform tablosundan yalnızca prob ID ve gen sembolünü, dosyayı belleğe almadan (zcat + awk) çıkarır.
# Clariom D gibi büyük platformların tam tablosu birkaç GB yer kaplayabildiği için bu şekilde okunur.
platform_sembol_oku <- function(dosya) {
  baslik <- system(sprintf("zcat '%s' | grep -m1 -A1 '^!platform_table_begin' | tail -n 1", dosya), intern = TRUE)
  if (!length(baslik)) return(NULL)
  sut <- strsplit(baslik, "\t")[[1]]
  aday <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "gene_symbol", "GeneSymbol")
  if (any(aday %in% sut)) { i <- match(aday[aday %in% sut][1], sut); tur <- "sembol"
  } else if ("gene_assignment" %in% sut) { i <- match("gene_assignment", sut); tur <- "ga"
  } else return(NULL)
  awk <- if (tur == "ga")
    sprintf("f&&/^!platform_table_end/{exit} f{n=split($%d,a,\" // \"); print $1 \"\\t\" (n>=2?a[2]:\"\")} /^!platform_table_begin/{f=1; getline}", i)
  else
    sprintf("f&&/^!platform_table_end/{exit} f{s=$%d; sub(/ \\/\\/\\/.*/,\"\",s); print $1 \"\\t\" s} /^!platform_table_begin/{f=1; getline}", i)
  komut <- sprintf("zcat '%s' | awk -F'\\t' '%s'", dosya, awk)
  tab <- read.delim(pipe(komut), header = FALSE, quote = "", comment.char = "",
                    col.names = c("ID", "sembol"), colClasses = "character")
  tab$sembol <- trimws(tab$sembol)
  tab
}

# GEOquery'nin yeni sürümleri (Bioconductor 3.24+) SummarizedExperiment döndürür; eski sürümler
# ExpressionSet döndürür. İki durumda da ExpressionSet'e çevrilir.
geo_oku <- function(f) {
  es <- tryCatch(getGEO(filename = f, getGPL = FALSE, returnType = "ExpressionSet"),
                 error = function(e) getGEO(filename = f, getGPL = FALSE))
  if (inherits(es, "SummarizedExperiment")) {
    m <- as.matrix(SummarizedExperiment::assay(es))
    pd <- as.data.frame(SummarizedExperiment::colData(es), optional = TRUE)
    rownames(pd) <- colnames(m)
    es <- ExpressionSet(assayData = m, phenoData = AnnotatedDataFrame(pd))
  }
  pd <- pData(es)
  if (!any(grepl(":ch1$", colnames(pd)))) {        # "anahtar: değer" sütunlarından :ch1 sütunları üret
    for (c in grep("^characteristics", colnames(pd), value = TRUE)) {
      v <- as.character(pd[[c]]); k <- unique(sub(":.*$", "", v[grepl(":", v)]))
      if (length(k) == 1) pd[[paste0(k, ":ch1")]] <- trimws(sub("^[^:]*:", "", v))
    }
    pData(es) <- pd
  }
  es
}

gse_yukle <- function(acc) {
  message("\n>>> ", acc, " indiriliyor...")
  url <- paste0(ftp_klasor(acc, "series"), "matrix/")
  yerel <- Filter(Negate(is.null), list(indir(paste0(url, acc, "_series_matrix.txt.gz"))))
  if (!length(yerel)) {   # birden çok platformlu seriler: dosya adları -GPLxxx içerir
    dosyalar <- suppressWarnings(dizin_listesi(url, paste0(acc, "[^\"<>/]*_series_matrix\\.txt\\.gz")))
    yerel <- Filter(Negate(is.null), lapply(paste0(url, dosyalar), indir))
  }
  if (!length(yerel)) stop(acc, ": seri matrisi indirilemedi (bağlantı tanı mesajlarına bakın).")
  liste <- lapply(yerel, geo_oku)
  en_buyuk <- which.max(sapply(liste, ncol))
  es <- liste[[en_buyuk]]
  if (nrow(exprs(es)) == 0)
    stop(acc, ": seri matrisinde ifade verisi yok (muhtemelen RNA-seq). Bu pilot yalnızca mikroarray içindir.")

  bas <- readLines(gzfile(yerel[[en_buyuk]]), n = 200, warn = FALSE)
  gpl_id <- gsub("[^A-Z0-9]", "", sub("^!Series_platform_id\t", "", grep("^!Series_platform_id", bas, value = TRUE)[1]))
  message("    Platform: ", gpl_id, " | ", nrow(es), " prob x ", ncol(es), " örnek")
  gpl <- NULL
  for (u in c(paste0(ftp_klasor(gpl_id, "platforms"), "annot/", gpl_id, ".annot.gz"),
              paste0(ftp_klasor(gpl_id, "platforms"), "soft/", gpl_id, "_family.soft.gz"))) {
    f <- indir(u)
    if (!is.null(f)) gpl <- platform_sembol_oku(f)
    if (!is.null(gpl) && mean(gpl$sembol != "" & gpl$sembol != "---") > 0.1) break
  }
  if (is.null(gpl)) stop(acc, ": platform anotasyonu (", gpl_id, ") indirilemedi veya gen sembolü sütunu yok.")
  s <- gpl$sembol[match(rownames(es), gpl$ID)]
  rm(gpl); invisible(gc())
  message("    Gen sembolü eşlenen prob oranı: ", round(100 * mean(!is.na(s) & s != "" & s != "---"), 1), "%")

  x <- exprs(es)
  if (max(x, na.rm = TRUE) > 100) x <- log2(pmax(x, 0) + 1)   # log2 değilse dönüştür
  ok <- !is.na(s) & s != "" & s != "---"
  x <- x[ok, , drop = FALSE]; s <- s[ok]
  # Aynı gene ait birden çok prob: ortalaması en yüksek prob seçilir
  if (identical(PROB_BIRLESTIRME, "ortalama")) {
    x <- rowsum(x, s, na.rm = TRUE) / as.vector(table(s)[sort(unique(s))])
  } else {
    ort <- rowMeans(x, na.rm = TRUE)
    sira <- order(s, -ort)
    x <- x[sira, , drop = FALSE]; s <- s[sira]
    x <- x[!duplicated(s), , drop = FALSE]; rownames(x) <- s[!duplicated(s)]
  }
  message("    Analize giren gen sayısı: ", nrow(x), " (prob birleştirme: ", PROB_BIRLESTIRME, ")")
  list(acc = acc, x = x, pd = pData(es))
}

# RNA-seq verisi: ifade NCBI'nin GEO için hazırladığı ham sayım matrisinden, örnek bilgileri seri matrisinden okunur.
rnaseq_yukle <- function(acc) {
  message("\n>>> ", acc, " (RNA-seq) indiriliyor...")
  f <- indir(paste0(ftp_klasor(acc, "series"), "matrix/", acc, "_series_matrix.txt.gz"))
  if (is.null(f)) stop(acc, ": seri matrisi (örnek bilgileri) indirilemedi.")
  pd <- pData(geo_oku(f))
  kok <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
  say_ad <- paste0(acc, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"); ann_ad <- "Human.GRCh38.p13.annot.tsv.gz"
  say_f <- file.path(INDIRME_KLASORU, say_ad); ann_f <- file.path(INDIRME_KLASORU, ann_ad)
  gecerli <- function(dosya, ilk) file.exists(dosya) &&
    isTRUE(tryCatch(startsWith(readLines(gzfile(dosya), n = 1), ilk), error = function(e) FALSE))
  indir_dene <- function(url, hedef) tryCatch(download.file(url, hedef, mode = "wb", quiet = TRUE),
    error = function(e) message("    indirilemedi: ", url, " (", conditionMessage(e), ")"),
    warning = function(w) message("    indirilemedi: ", url, " (", conditionMessage(w), ")"))
  if (!gecerli(say_f, "GeneID")) indir_dene(paste0(kok, "&acc=", acc, "&format=file&file=", say_ad), say_f)
  if (!gecerli(say_f, "GeneID"))
    stop(acc, ": NCBI sayım dosyası alınamadı (", say_ad, "). Sunucu erişimi engellemiş olabilir.")
  for (u in c(paste0(kok, "&acc=", acc, "&format=file&file=", ann_ad), paste0(kok, "&format=file&file=", ann_ad)))
    if (!gecerli(ann_f, "GeneID")) indir_dene(u, ann_f)
  say <- read.delim(gzfile(say_f), check.names = FALSE, row.names = 1)
  if (gecerli(ann_f, "GeneID")) {
    ann <- read.delim(gzfile(ann_f), quote = "", check.names = FALSE, comment.char = "")
    sembol <- as.character(ann$Symbol[match(rownames(say), ann$GeneID)])
  } else {                                     # yedek: Bioconductor org.Hs.eg.db ile Entrez -> sembol
    message("    NCBI anotasyonu alınamadı; org.Hs.eg.db ile eşleniyor (ilk seferde kurulum birkaç dakika sürer)...")
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
    sembol <- unname(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(say),
                                                        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")))
  }
  say <- as.matrix(say); ortak <- intersect(colnames(say), rownames(pd))
  if (length(ortak) < 4) stop(acc, ": sayım matrisi ile örnek bilgileri eşleşmedi.")
  say <- say[, ortak]; pd <- pd[ortak, , drop = FALSE]
  kutup <- colSums(say)
  if (exists("RNASEQ_TMM") && isTRUE(RNASEQ_TMM)) {   # duyarlılık: edgeR TMM normalizasyonu
    if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR", update = FALSE, ask = FALSE)
    kutup <- kutup * edgeR::calcNormFactors(edgeR::DGEList(say))$samples$norm.factors
    message("    TMM normalizasyonu uygulandı")
  }
  cpm <- t(t(say) / kutup * 1e6)
  tut <- rowMeans(cpm > 1) >= 0.2 & !is.na(sembol) & sembol != ""
  x <- log2(cpm[tut, , drop = FALSE] + 1); s <- sembol[tut]
  ort <- rowMeans(x); sira <- order(s, -ort); x <- x[sira, , drop = FALSE]; s <- s[sira]
  x <- x[!duplicated(s), , drop = FALSE]; rownames(x) <- s[!duplicated(s)]
  message("    ", ncol(x), " örnek | analize giren gen sayısı: ", nrow(x), " (log2 CPM)")
  list(acc = acc, x = x, pd = pd, rnaseq = TRUE)
}

## ---------------------------------------------------------------------------
## 4) METADATA: ZAMAN VE KİŞİ SÜTUNLARINI BULMA
## ---------------------------------------------------------------------------
metadata_ozeti_yaz <- function(acc, pd) {
  sut <- grep(":ch1$|^characteristics|^title$|^source_name", colnames(pd), value = TRUE)
  dosya <- yol("01_metadata_", acc, ".txt")
  sink(dosya)
  cat("Veri seti:", acc, " | örnek sayısı:", nrow(pd), "\n\n")
  for (c in sut) {
    u <- unique(as.character(pd[[c]]))
    cat("SÜTUN:", c, " (", length(u), " farklı değer)\n  ", paste(head(u, 8), collapse = " | "),
        if (length(u) > 8) " | ..." else "", "\n\n")
  }
  sink()
  dosya
}

sutunlari_algila <- function(acc, pd) {
  if (!is.null(ELLE_AYAR[[acc]])) return(ELLE_AYAR[[acc]])
  sut <- grep(":ch1$", colnames(pd), value = TRUE)
  zaman_deseni <- "before|after|baseline|pre|post|year|month|week|visit|time|t0|t1|t2|t3|follow|maint|cid|surgery|control|lean|never|non.?obese"
  zaman <- NULL; kisi <- NULL
  for (c in sut) {
    v <- tolower(as.character(pd[[c]])); n <- length(unique(v))
    if (is.null(zaman) && n >= 2 && n <= 5 && mean(grepl(zaman_deseni, v)) > 0.6) zaman <- c
  }
  for (c in sut) {
    n <- length(unique(pd[[c]]))
    if (is.null(zaman) && n >= 2 && n <= 5 && grepl("time|visit|point|stage|cid|state|status|group|condition", tolower(c))) zaman <- c
  }
  for (c in sut) {
    v <- as.character(pd[[c]]); n <- length(unique(v))
    if (is.null(kisi) && grepl("subject|patient|individual|donor|participant|^id|code|pair", tolower(c)) &&
        n < length(v) && n > 2)
      kisi <- c
  }
  if (is.null(kisi)) {   # başlıktan kişi kimliği çıkarmayı dene (ör. "Subject 12, before")
    t <- as.character(pd$title)
    k <- sub(".*?(subject|patient|individual|id|p)[ _#-]*([0-9]+).*", "\\2", t, ignore.case = TRUE, perl = TRUE)
    if (mean(k != t) > 0.8 && length(unique(k)) < length(k)) { pd$..kisi <- k; kisi <- "..kisi" }
  }
  if (is.null(zaman)) return(NULL)
  if (is.null(kisi) && !is.null(pd$title)) {   # başlıktan zaman ifadelerini silip kalanı kişi kimliği say
    t <- tolower(as.character(pd$title))
    for (z in unique(tolower(as.character(pd[[zaman]])))) t <- gsub(z, " ", t, fixed = TRUE)
    t <- gsub("\\b(t[0-9]+|m[0-9]+|y[0-9]+|v[0-9]+|before|after|pre|post|baseline|follow.?up|surgery|op|month[s]?|year[s]?|[0-9]+ ?(m|mo|y|yr))\\b", " ", t, perl = TRUE)
    t <- gsub("[^a-z0-9]+", "_", t); t <- gsub("^_+|_+$", "", t)
    tab <- table(t)
    if (length(tab) < length(t) && length(tab) >= 3 && mean(tab >= 2) >= 0.5) { pd$..kisi <- t; kisi <- "..kisi" }
  }
  eslesmesiz <- is.null(kisi)
  if (eslesmesiz) { pd$..kisi <- rownames(pd); kisi <- "..kisi" }   # eşleştirme yok: her örnek ayrı kişi
  v <- as.character(pd[[zaman]])
  kon <- unique(v[grepl("control|lean|never|non.?obese", v, ignore.case = TRUE)])
  aday <- setdiff(unique(v), kon)
  bas <- aday[grepl("before|baseline|pre|start|^t0|t0$|cid ?1|week 0|^0 ?(month|year)", aday, ignore.case = TRUE)]
  list(zaman = zaman, kisi = kisi, baslangic = if (length(bas)) bas[1] else NULL,
       kontrol = if (length(kon) == 1) kon else NULL, eslesmesiz = eslesmesiz, pd = pd)
}

## ---------------------------------------------------------------------------
## 5) SKORLAR
## ---------------------------------------------------------------------------
# ssGSEA (Barbie ve ark. 2009; GSVA paketindeki uygulamayla birebir aynı sonuç verir)
ssgsea_skor <- function(x, setler, alpha = 0.25) {
  x <- x[rowSums(is.na(x)) == 0, , drop = FALSE]
  setler <- lapply(setler, function(g) intersect(g, rownames(x)))
  setler <- setler[sapply(setler, length) >= 5]
  sk <- sapply(seq_len(ncol(x)), function(j) {
    r <- rank(x[, j], ties.method = "average")
    o <- order(r, decreasing = TRUE); ro <- r[o]; gen_o <- rownames(x)[o]
    sapply(setler, function(g) {
      hit <- gen_o %in% g
      w <- ro^alpha * hit
      sum(cumsum(w) / sum(w) - cumsum(!hit) / sum(!hit))
    })
  })
  sk <- matrix(sk, nrow = length(setler), dimnames = list(names(setler), colnames(x)))
  sk / diff(range(sk))
}
z_skor <- function(x, genler) {          # basit, yorumlanması kolay ikinci skor
  g <- intersect(genler, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

## ---------------------------------------------------------------------------
## 5b) SenMayo GENLERİNİN SABİT HÜCRE TİPİ SINIFLANDIRMASI (Human Protein Atlas)
## ---------------------------------------------------------------------------
# Her SenMayo geni, HPA tek hücre RNA verisinde bağışıklık hücrelerinde mi yoksa yağ dokusunun
# stromal/adiposit hücrelerinde mi daha çok ifade edildiğine göre BİR KEZ sınıflandırılır ve aynı
# liste tüm kohortlarda kullanılır.
HPA_URLLERI <- c("https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
                 "https://v24.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
                 "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip",
                 "https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip")
IMMUN_DESEN   <- "macrophage|monocyte|t-cell|t cell|b-cell|b cell|nk-cell|nk cell|mast|dendritic|granulocyte|neutrophil|plasma|kupffer|langerhans|hofbauer|microglial"
STROMAL_DESEN <- "adipocyte|fibroblast|endothelial|smooth muscle|mesothelial|pericyte"

hpa_siniflandir <- function(genler) {
  z <- file.path(INDIRME_KLASORU, "hpa_rna_single_cell_type.tsv.zip")
  if (!file.exists(z) || file.size(z) < 2000) {
    for (u in HPA_URLLERI) {
      ok <- tryCatch({ download.file(u, z, mode = "wb", quiet = TRUE); file.size(z) > 2000 },
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) { message("    HPA tek hücre verisi indirildi: ", u); break }
      message("    HPA indirilemedi: ", u)
    }
  }
  if (!file.exists(z) || file.size(z) < 2000) return(NULL)
  icerik <- unzip(z, list = TRUE)$Name[1]
  h <- read.delim(unz(z, icerik), check.names = FALSE, stringsAsFactors = FALSE)
  cn <- tolower(colnames(h))
  gad <- colnames(h)[match("gene name", cn)]
  if (is.na(gad)) { message("    HPA dosyasında 'Gene name' sütunu yok; sınıflandırma atlandı."); return(NULL) }
  h <- h[h[[gad]] %in% genler, , drop = FALSE]
  if ("cell type" %in% cn) {                                   # uzun format
    ct <- colnames(h)[match("cell type", cn)]; nt <- colnames(h)[match("ntpm", cn)]
    m <- tapply(as.numeric(h[[nt]]), list(h[[gad]], h[[ct]]), max)
  } else {                                                     # geniş format
    say <- sapply(h, is.numeric); m <- as.matrix(h[, say, drop = FALSE]); rownames(m) <- h[[gad]]
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
  }
  m[is.na(m)] <- 0
  ct <- tolower(colnames(m)); imm <- grepl(IMMUN_DESEN, ct); str <- grepl(STROMAL_DESEN, ct)
  if (sum(imm) < 3 || sum(str) < 2) { message("    HPA hücre tipleri tanınamadı; sınıflandırma atlandı."); return(NULL) }
  imm_max <- apply(m[, imm, drop = FALSE], 1, max); str_max <- apply(m[, str, drop = FALSE], 1, max)
  pay <- imm_max / (imm_max + str_max + 1e-9)
  data.frame(gen = rownames(m), immun_max_nTPM = round(imm_max, 1), stromal_max_nTPM = round(str_max, 1),
             immun_payi = round(pay, 3),
             sinif = ifelse(imm_max + str_max < 1, "dusuk_ifade",
                     ifelse(pay >= 0.6, "immun", ifelse(pay <= 0.4, "stromal_adiposit", "karisik"))),
             stringsAsFactors = FALSE)
}
SABIT_SINIF <- NULL   # çalıştırma bölümünde doldurulur

## ---------------------------------------------------------------------------
## 6) ANALİZ (bir veri seti için)
## ---------------------------------------------------------------------------
analiz_et <- function(d, setler) {
  acc <- d$acc
  metadata_ozeti_yaz(acc, d$pd)
  a <- sutunlari_algila(acc, d$pd)
  if (is.null(a) || is.null(a$baslangic)) {
    message("!!! ", acc, ": zaman/kişi sütunları otomatik bulunamadı. ",
            yol("01_metadata_", acc, ".txt"), " dosyasına bakıp ELLE_AYAR'ı doldurun.")
    return(NULL)
  }
  pd <- if (!is.null(a$pd)) a$pd else d$pd
  sutun_bul <- function(ad) {                      # sütun adı sürümler arasında küçük farklılık gösterebilir
    if (is.null(ad) || ad %in% colnames(pd)) return(ad)
    aday <- c(sub(":ch1$", "", ad), make.names(ad), make.names(sub(":ch1$", "", ad)))
    b <- aday[aday %in% colnames(pd)]
    if (length(b)) return(b[1])
    stop("'", ad, "' sütunu bulunamadı. Mevcut sütunlar: ", paste(colnames(pd), collapse = ", "))
  }
  a$zaman <- sutun_bul(a$zaman); a$kisi <- sutun_bul(a$kisi); a$grup <- sutun_bul(a$grup)
  zaman <- factor(as.character(pd[[a$zaman]]))
  dz <- setdiff(levels(zaman), c(a$baslangic, a$kontrol))            # sıra: başlangıç, takip(ler), kontrol
  zaman <- factor(zaman, levels = c(a$baslangic, dz, a$kontrol))
  kisi  <- factor(as.character(pd[[a$kisi]]))
  grup  <- if (!is.null(a$grup)) as.character(pd[[a$grup]]) else rep(NA_character_, nrow(pd))
  message("    Zaman sütunu: ", a$zaman, " (", paste(levels(zaman), collapse = ", "), ")")
  message("    Kişi sütunu : ", a$kisi, " (", nlevels(kisi), " kişi)")

  # --- Makrofaj skoru ve SenMayo'nun bağışıklık/bağışıklık dışı bölünmesi ---
  # Makrofaj seti SenMayo ile hiç ortak gen içermez (döngüsellik olmasın diye).
  setler$Makrofaj <- setdiff(MAKROFAJ_GENLERI, setler$SenMayo)
  setler$T_hucresi <- setdiff(T_HUCRE_GENLERI, setler$SenMayo)
  setler$Mast_hucresi <- setdiff(MAST_GENLERI, setler$SenMayo)
  setler$Adiposit <- setdiff(ADIPOSIT_GENLERI, setler$SenMayo)
  g_sen_tum <- intersect(setler$SenMayo, rownames(d$x))
  if (!is.null(SABIT_SINIF)) {                    # HPA'ya dayalı sabit sınıflandırma (tüm kohortlarda aynı)
    setler$SenMayo_immun      <- intersect(SABIT_SINIF$gen[SABIT_SINIF$sinif == "immun"], g_sen_tum)
    setler$SenMayo_immun_disi <- intersect(SABIT_SINIF$gen[SABIT_SINIF$sinif == "stromal_adiposit"], g_sen_tum)
    kaynak <- "HPA tek hücre verisi"
  } else {                                        # yedek: ilk kohortta korelasyonla belirlenip sabitlenir
    imm <- ssgsea_skor(d$x, setler[c("Makrofaj", "T_hucresi", "Mast_hucresi")])
    r_imm <- sapply(rownames(imm), function(k)
      apply(d$x[g_sen_tum, , drop = FALSE], 1, function(v) suppressWarnings(cor(v, imm[k, ]))))
    r_imm <- matrix(r_imm, nrow = length(g_sen_tum), dimnames = list(g_sen_tum, rownames(imm)))
    r_max <- apply(r_imm, 1, function(r) if (all(is.na(r))) NA else max(r, na.rm = TRUE))
    immun_mu <- !is.na(r_max) & r_max >= MAKROFAJ_KOR_ESIGI
    SABIT_SINIF <<- data.frame(gen = g_sen_tum, sinif = ifelse(immun_mu, "immun", "stromal_adiposit"))
    write.csv(SABIT_SINIF, yol("06_SenMayo_siniflandirma.csv"), row.names = FALSE)
    setler$SenMayo_immun <- g_sen_tum[immun_mu]; setler$SenMayo_immun_disi <- g_sen_tum[!immun_mu]
    kaynak <- paste0(acc, " korelasyonu (yedek yöntem; sonraki kohortlarda aynı liste kullanılır)")
  }
  message("    SenMayo sınıflandırması (", kaynak, "): ", length(setler$SenMayo_immun), " bağışıklık, ",
          length(setler$SenMayo_immun_disi), " stromal/adiposit geni")

  kapsam <- sapply(setler, function(g) sum(g %in% rownames(d$x)))
  sk <- ssgsea_skor(d$x, setler)
  df <- data.frame(ornek = colnames(d$x), kisi = kisi, zaman = zaman, grup = grup,
                   t(sk), SenMayo_z = z_skor(d$x, setler$SenMayo), check.names = FALSE)
  # Makrofaj içeriğinden arındırılmış SenMayo skoru (doğrusal regresyon artıkları)
  df$SenMayo_makrofaj_duzeltilmis <- unname(residuals(lm(df$SenMayo ~ df$Makrofaj)))
  # Klinik değişken (BMI) varsa ekle
  bmi_sut <- grep("(^|[^a-z])bmi([^a-z]|$)|body mass", colnames(pd), ignore.case = TRUE, value = TRUE)
  bmi_sut <- bmi_sut[!grepl("^characteristics", bmi_sut)]
  if (length(bmi_sut)) {
    deger <- as.character(pd[[bmi_sut[1]]])
    df$BMI <- suppressWarnings(as.numeric(sub(",", ".", gsub("[^0-9.,]", "", deger))))
    message("    BMI sütunu: ", bmi_sut[1], " (", sum(!is.na(df$BMI)), " dolu değer)")
  } else message("    BMI sütunu bulunamadı.")
  write.csv(df, yol("02_skorlar_", acc, ".csv"), row.names = FALSE)

  # --- Eşleştirilmiş karşılaştırmalar: her zaman noktası vs başlangıç ---
  skor_sut <- intersect(c(names(setler), "SenMayo_z", "SenMayo_makrofaj_duzeltilmis"), colnames(df))
  eslestir <- function(dd, grup_etiketi) {
    out <- list()
    for (zn in dz) {
      ortak <- intersect(dd$kisi[dd$zaman == a$baslangic], dd$kisi[dd$zaman == zn])
      for (s in skor_sut) {
        if (length(ortak) >= 3) {                       # eşleştirilmiş test
          once  <- dd[dd$zaman == a$baslangic, ]; sonra <- dd[dd$zaman == zn, ]
          once  <- once[match(ortak, once$kisi), s]; sonra <- sonra[match(ortak, sonra$kisi), s]
          fark <- sonra - once
          out[[length(out) + 1]] <- data.frame(
            veri_seti = acc, grup = grup_etiketi, karsilastirma = paste(zn, "vs", a$baslangic), skor = s,
            test = "eslestirilmis", n_cift = length(ortak), ort_degisim = mean(fark), cohen_dz = mean(fark) / sd(fark),
            p_eslestirilmis_t = t.test(sonra, once, paired = TRUE)$p.value,
            p_wilcoxon = suppressWarnings(wilcox.test(sonra, once, paired = TRUE)$p.value))
        } else {                                        # eşleştirme yoksa: bağımsız gruplar (Welch)
          x0 <- dd[dd$zaman == a$baslangic, s]; x1 <- dd[dd$zaman == zn, s]
          if (length(x0) < 3 || length(x1) < 3) next
          out[[length(out) + 1]] <- data.frame(
            veri_seti = acc, grup = grup_etiketi, karsilastirma = paste(zn, "vs", a$baslangic), skor = s,
            test = "eslestirilmemis", n_cift = min(length(x0), length(x1)), ort_degisim = mean(x1) - mean(x0),
            cohen_dz = (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2),
            p_eslestirilmis_t = t.test(x1, x0)$p.value,
            p_wilcoxon = suppressWarnings(wilcox.test(x1, x0)$p.value))
        }
      }
    }
    do.call(rbind, out)
  }
  sonuc <- eslestir(df, "tum_orneklem")
  if (!is.null(a$grup))                       # alt grup analizi (ör. diyet tipi)
    for (g in unique(na.omit(df$grup))) sonuc <- rbind(sonuc, eslestir(df[df$grup %in% g, ], g))

  # --- Hiç obez olmamış kontrollerle karşılaştırma (eşleştirilmemiş, Welch t) ---
  kontrol_sonuc <- NULL
  if (!is.null(a$kontrol) && sum(df$zaman == a$kontrol) >= 3) {
    for (zn in c(a$baslangic, dz)) for (s in skor_sut) {
      x1 <- df[df$zaman == zn, s]; x0 <- df[df$zaman == a$kontrol, s]
      kontrol_sonuc <- rbind(kontrol_sonuc, data.frame(
        veri_seti = acc, karsilastirma = paste(zn, "vs", a$kontrol), skor = s,
        n_grup = length(x1), n_kontrol = length(x0), ort_fark = mean(x1) - mean(x0),
        cohen_d = (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2),
        p_welch = t.test(x1, x0)$p.value))
    }
  }

  # --- Özgüllük: aynı büyüklükte rastgele gen setleri ---
  g_sen <- intersect(setler$SenMayo, rownames(d$x))
  ozgulluk <- list()
  for (zn in dz) {
    ortak <- intersect(df$kisi[df$zaman == a$baslangic], df$kisi[df$zaman == zn])
    if (length(ortak) < 3) next
    i0 <- match(ortak, df$kisi[df$zaman == a$baslangic]); i1 <- match(ortak, df$kisi[df$zaman == zn])
    s0 <- which(df$zaman == a$baslangic)[i0];             s1 <- which(df$zaman == zn)[i1]
    degisim <- function(genler) { z <- z_skor(d$x, genler); mean(z[s1] - z[s0]) }
    gercek <- degisim(g_sen)
    rast <- replicate(RASTGELE_SET_SAYISI, degisim(sample(rownames(d$x), length(g_sen))))
    ozgulluk[[zn]] <- data.frame(veri_seti = acc, karsilastirma = paste(zn, "vs", a$baslangic),
      gercek_degisim = gercek, rastgele_ortalama = mean(rast),
      ampirik_p_iki_yonlu = (sum(abs(rast - mean(rast)) >= abs(gercek - mean(rast))) + 1) / (RASTGELE_SET_SAYISI + 1))
  }
  ozgulluk <- do.call(rbind, ozgulluk)

  # --- Gen düzeyi: SenMayo genleri (limma, kişi bloklu eşleştirilmiş tasarım) ---
  tutulan <- kisi %in% names(which(table(kisi) >= 2))
  if (sum(tutulan) >= 6) {
    K <- droplevels(kisi[tutulan]); Z <- droplevels(zaman[tutulan])
    tasarim <- model.matrix(~ K + Z)
  } else {                                   # eşleştirme yok: kontrol grubu hariç basit grup karşılaştırması
    tutulan <- !is.na(zaman) & !(zaman %in% a$kontrol)
    Z <- droplevels(zaman[tutulan]); tasarim <- model.matrix(~ Z)
  }
  fit <- eBayes(lmFit(d$x[, tutulan], tasarim), trend = isTRUE(d$rnaseq))
  zaman_kat <- grep("^Z", colnames(tasarim))
  genler <- do.call(rbind, lapply(zaman_kat, function(k) {
    tt <- topTable(fit, coef = k, number = Inf)
    tt$gen <- rownames(tt)
    tt$karsilastirma <- paste(sub("^Z", "", colnames(tasarim)[k]), "vs", a$baslangic)
    tt$SenMayo_geni <- tt$gen %in% g_sen
    tt
  }))
  write.csv(genler[genler$SenMayo_geni, ], yol("04_SenMayo_genleri_", acc, ".csv"), row.names = FALSE)
  temel <- genler[genler$gen %in% TEMEL_SENESANS_GENLERI, c("gen", "karsilastirma", "logFC", "P.Value", "adj.P.Val")]
  write.csv(temel[order(temel$karsilastirma, temel$gen), ], yol("04_temel_senesans_genleri_", acc, ".csv"), row.names = FALSE)

  # --- Grafik ---
  p <- ggplot(df, aes(zaman, SenMayo, group = kisi)) +
    geom_line(alpha = 0.25, colour = "grey40") +
    geom_boxplot(aes(group = zaman), width = 0.35, fill = "#7FA7C9", alpha = 0.6, outlier.shape = NA) +
    geom_point(size = 1, alpha = 0.6) +
    labs(title = paste0(acc, ": SenMayo ssGSEA skoru"), x = NULL, y = "SenMayo skoru") +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  if (!is.null(a$grup)) p <- p + facet_wrap(~ grup)
  ggsave(yol("05_grafik_", acc, ".png"), p, width = if (is.null(a$grup)) 6.5 else 9, height = 4.8, dpi = 300)

  # --- Doku bileşimine göre düzeltilmiş değişim ---
  # Eşleştirilmiş farklar üzerinde: Δskor = b0 + b1·ΔAdiposit + b2·ΔMakrofaj. b0, doku bileşimi
  # (adiposit ve makrofaj payı) değişmeseydi beklenecek skor değişimidir.
  bilesim <- NULL
  for (zn in dz) {
    ortak <- intersect(df$kisi[df$zaman == a$baslangic], df$kisi[df$zaman == zn])
    if (length(ortak) < 10 || !all(c("Adiposit", "Makrofaj") %in% colnames(df))) next
    o <- df[df$zaman == a$baslangic, ][match(ortak, df$kisi[df$zaman == a$baslangic]), ]
    s2 <- df[df$zaman == zn, ][match(ortak, df$kisi[df$zaman == zn]), ]
    for (sk_ad in intersect(c("SenMayo", "SenMayo_immun", "SenMayo_immun_disi", "Fridman_SEN_UP"), colnames(df))) {
      dd <- data.frame(y = s2[[sk_ad]] - o[[sk_ad]], adip = s2$Adiposit - o$Adiposit, mak = s2$Makrofaj - o$Makrofaj)
      m0 <- summary(lm(y ~ 1, dd))$coefficients; m1 <- summary(lm(y ~ adip + mak, dd))$coefficients
      bilesim <- rbind(bilesim, data.frame(veri_seti = acc, karsilastirma = paste(zn, "vs", a$baslangic),
        skor = sk_ad, n_cift = nrow(dd), ham_degisim = m0[1, 1], ham_se = m0[1, 2], ham_p = m0[1, 4],
        duzeltilmis_degisim = m1[1, 1], duzeltilmis_se = m1[1, 2], duzeltilmis_p = m1[1, 4],
        adiposit_katsayi_p = m1["adip", 4], makrofaj_katsayi_p = m1["mak", 4]))
    }
  }

  # --- Duyarlılık analizleri: (1) farklı HPA eşikleri, (2) tek gen çıkarma ---
  # Not: ssGSEA'nın aralık normalizasyonu Cohen dz'yi etkilemediği için burada atlanır.
  duy <- NULL
  if (isTRUE(DUYARLILIK)) {
    message("    Duyarlılık analizleri hesaplanıyor...")
    gx <- d$x
    sir <- lapply(seq_len(ncol(gx)), function(j) {
      r <- rank(gx[, j], ties.method = "average"); o <- order(r, decreasing = TRUE)
      list(g = rownames(gx)[o], w = r[o]^0.25) })
    es <- function(gs) vapply(sir, function(sv) {
      hit <- sv$g %in% gs; w <- sv$w * hit
      sum(cumsum(w) / sum(w) - cumsum(!hit) / sum(!hit)) }, 0)
    esl_dz <- function(sc, zn) {
      ortak <- intersect(df$kisi[df$zaman == a$baslangic], df$kisi[df$zaman == zn])
      if (length(ortak) < 3) return(c(NA, NA, NA))
      o <- sc[df$zaman == a$baslangic][match(ortak, df$kisi[df$zaman == a$baslangic])]
      s3 <- sc[df$zaman == zn][match(ortak, df$kisi[df$zaman == zn])]
      f <- s3 - o; c(length(ortak), mean(f) / sd(f), t.test(f)$p.value) }
    satir <- function(analiz, bil, ng, zn, r, rr = NULL) data.frame(
      veri_seti = acc, analiz = analiz, bilesen = bil, n_gen = ng, karsilastirma = paste(zn, "vs", a$baslangic),
      n_cift = r[1], cohen_dz = r[2], p = r[3],
      dz_min = if (is.null(rr)) NA else min(rr[, 2]), dz_max = if (is.null(rr)) NA else max(rr[, 2]),
      p_en_buyuk = if (is.null(rr)) NA else max(rr[, 3]))
    if (!is.null(SABIT_SINIF$immun_payi)) {
      esikler <- list(ana = c(0.6, 0.4), siki = c(0.7, 0.3), gevsek = c(0.5, 0.5))
      ok <- SABIT_SINIF$sinif != "dusuk_ifade"; pay <- SABIT_SINIF$immun_payi; gn <- SABIT_SINIF$gen
      for (nm in names(esikler)) {
        e <- esikler[[nm]]
        kume <- list(immun = intersect(gn[ok & pay >= e[1]], rownames(gx)),
                     stromal = intersect(gn[ok & (if (nm == "gevsek") pay < e[2] else pay <= e[2])], rownames(gx)))
        for (bil in names(kume)) {
          if (length(kume[[bil]]) < 5) next
          sc <- es(kume[[bil]])
          for (zn in dz) { r <- esl_dz(sc, zn); if (!is.na(r[1])) duy <- rbind(duy, satir(paste0("esik_", nm), bil, length(kume[[bil]]), zn, r)) }
        }
      }
    }
    for (bil in c("immun", "stromal")) {
      gs <- if (bil == "immun") setler$SenMayo_immun else setler$SenMayo_immun_disi
      gs <- intersect(gs, rownames(gx)); if (length(gs) < 6) next
      tam <- es(gs); loo <- lapply(gs, function(g) es(setdiff(gs, g)))
      for (zn in dz) {
        r <- esl_dz(tam, zn); if (is.na(r[1])) next
        rr <- t(sapply(loo, esl_dz, zn = zn))
        duy <- rbind(duy, satir("tek_gen_cikarma", bil, length(gs), zn, r, rr))
      }
    }
  }

  # BMI değişimi ile skor değişimi ilişkisi (varsa)
  bmi_kor <- NULL
  if (!is.null(df$BMI) && sum(!is.na(df$BMI)) > 10) {
    for (zn in dz) {
      ortak <- intersect(df$kisi[df$zaman == a$baslangic], df$kisi[df$zaman == zn])
      if (length(ortak) < 8) next
      o <- df[df$zaman == a$baslangic, ][match(ortak, df$kisi[df$zaman == a$baslangic]), ]
      s <- df[df$zaman == zn, ][match(ortak, df$kisi[df$zaman == zn]), ]
      dbmi <- 100 * (s$BMI - o$BMI) / o$BMI
      for (sk_ad in c("SenMayo", "SenMayo_immun", "SenMayo_immun_disi", "Makrofaj")) {
        k <- suppressWarnings(cor.test(dbmi, s[[sk_ad]] - o[[sk_ad]], method = "spearman"))
        bmi_kor <- rbind(bmi_kor, data.frame(veri_seti = acc, karsilastirma = paste(zn, "vs", a$baslangic),
                         baslangic_BMI_ort = mean(o$BMI, na.rm = TRUE), BMI_degisim_yuzde_ort = mean(dbmi, na.rm = TRUE),
                         skor = sk_ad, spearman_rho = unname(k$estimate), p = k$p.value))
      }
    }
  }
  list(sonuc = sonuc, ozgulluk = ozgulluk, kapsam = kapsam, bmi = bmi_kor, kontrol = kontrol_sonuc, bilesim = bilesim, duy = duy)
}

## ---------------------------------------------------------------------------
## 7) ÇALIŞTIR
## ---------------------------------------------------------------------------
setler <- gen_setlerini_al()
message("Gen setleri: ", paste(names(setler), sapply(setler, length), sep = "=", collapse = ", "))
message("\n>>> SenMayo genleri Human Protein Atlas tek hücre verisine göre sınıflandırılıyor...")
SABIT_SINIF <- tryCatch(hpa_siniflandir(setler$SenMayo), error = function(e) { message("    HPA hatası: ", conditionMessage(e)); NULL })
if (!is.null(SABIT_SINIF)) {
  write.csv(SABIT_SINIF, yol("06_SenMayo_siniflandirma.csv"), row.names = FALSE)
  message("    ", paste(names(table(SABIT_SINIF$sinif)), table(SABIT_SINIF$sinif), sep = "=", collapse = ", "))
} else message("    HPA kullanılamadı; yedek yöntem (ilk kohortta korelasyon) kullanılacak.")

tum <- list()
for (acc in VERI_SETLERI) {
  tum[[acc]] <- tryCatch(analiz_et(if (acc %in% RNASEQ_SETLERI) rnaseq_yukle(acc) else gse_yukle(acc), setler),
                         error = function(e) { message("!!! ", acc, " HATA: ", conditionMessage(e)); NULL })
}
tum <- tum[!sapply(tum, is.null)]

if (length(tum)) {
  sonuc    <- do.call(rbind, lapply(tum, `[[`, "sonuc"))
  ozgulluk <- do.call(rbind, lapply(tum, `[[`, "ozgulluk"))
  bmi      <- do.call(rbind, lapply(tum, `[[`, "bmi"))
  kontrol  <- do.call(rbind, lapply(tum, `[[`, "kontrol"))
  bilesim  <- do.call(rbind, lapply(tum, `[[`, "bilesim"))
  duyarlilik <- do.call(rbind, lapply(tum, `[[`, "duy"))
  if (!is.null(duyarlilik)) write.csv(duyarlilik, yol("08_duyarlilik.csv"), row.names = FALSE)
  if (!is.null(bilesim)) write.csv(bilesim, yol("03_bilesim_duzeltmesi.csv"), row.names = FALSE)
  if (!is.null(kontrol)) write.csv(kontrol, yol("03_kontrol_karsilastirma.csv"), row.names = FALSE)
  write.csv(sonuc,    yol("03_eslestirilmis_testler.csv"), row.names = FALSE)
  write.csv(ozgulluk, yol("03_ozgulluk_rastgele_setler.csv"), row.names = FALSE)
  if (!is.null(bmi)) write.csv(bmi, yol("03_BMI_korelasyon.csv"), row.names = FALSE)

  sink(yol("00_OZET.txt"))
  cat("PİLOT ANALİZ ÖZETİ -", format(Sys.time()), "\n\n")
  for (n in names(tum)) cat(n, "- gen seti kapsamı:",
                            paste(names(tum[[n]]$kapsam), tum[[n]]$kapsam, sep = "=", collapse = ", "), "\n")
  ana_skorlar <- c("SenMayo", "SenMayo_immun", "SenMayo_immun_disi", "Makrofaj", "T_hucresi", "Mast_hucresi", "Adiposit", "Fridman_SEN_UP")
  cat("\nEşleştirilmiş testler (tüm örneklem):\n")
  ss <- sonuc[sonuc$skor %in% ana_skorlar & sonuc$grup == "tum_orneklem", ]
  ss <- ss[order(ss$veri_seti, ss$karsilastirma, match(ss$skor, ana_skorlar)), ]
  print(ss[, c("veri_seti","karsilastirma","skor","test","n_cift","ort_degisim","cohen_dz","p_eslestirilmis_t")],
        row.names = FALSE, digits = 3)
  if (!is.null(kontrol)) {
    cat("\nHiç obez olmamış kontrollere göre:\n")
    kk <- kontrol[kontrol$skor %in% ana_skorlar, ]
    kk <- kk[order(kk$karsilastirma, match(kk$skor, ana_skorlar)), ]
    print(kk[, c("veri_seti","karsilastirma","skor","n_grup","n_kontrol","ort_fark","cohen_d","p_welch")],
          row.names = FALSE, digits = 3)
  }
  cat("\nÖzgüllük (rastgele gen setlerine göre):\n"); print(ozgulluk, row.names = FALSE, digits = 3)
  if (!is.null(bilesim)) {
    cat("\nDoku bileşimine (adiposit + makrofaj) göre düzeltilmiş değişim:\n")
    print(bilesim, row.names = FALSE, digits = 3)
  }
  if (!is.null(SABIT_SINIF)) { cat("\nSenMayo sınıflandırması:\n"); print(table(SABIT_SINIF$sinif)) }
  if (!is.null(duyarlilik)) {
    cat("\nDuyarlılık analizleri (Cohen dz; tek gen çıkarmada dz aralığı ve en büyük p):\n")
    print(duyarlilik, row.names = FALSE, digits = 3)
  }
  if (!is.null(bmi)) { cat("\nBMI değişimi (%) ile skor değişimi (Spearman):\n"); print(bmi, row.names = FALSE, digits = 3) }
  for (f in list.files(CIKTI_KLASORU, pattern = "^04_temel", full.names = TRUE)) {
    cat("\nTemel senesans genleri -", sub("04_temel_senesans_genleri_|\\.csv", "", basename(f)), ":\n")
    print(read.csv(f), row.names = FALSE, digits = 3)
  }
  cat("\nYorum: ort_degisim < 0 ve p < 0.05 ise senesans skoru kilo kaybından sonra düşmüştür.\n",
      "Özgüllük p < 0.05 ise düşüş rastgele gen setlerinden beklenenden büyüktür.\n",
      "Kontrol karşılaştırmasında takip noktasının p değeri büyükse skor sağlıklı düzeye yaklaşmıştır.\n",
      "Bileşim düzeltmesinde duzeltilmis_p < 0.05 ise değişim doku bileşimi değişiminden bağımsızdır.\n")
  sink()
  cat(readLines(yol("00_OZET.txt")), sep = "\n")
}

## ---------------------------------------------------------------------------
## 8) YAYIN ŞEKİLLERİ (07_sekil*.png / .pdf)
## ---------------------------------------------------------------------------
# Şekil etiketleri İngilizcedir (uluslararası dergi için)
KOHORT_ETIKET <- c(GSE65540 = "Surgery (GSE65540)", GSE59034 = "Surgery (GSE59034)",
                   GSE199063 = "Surgery (GSE199063)", GSE77962 = "Diet (GSE77962)")
SKOR_ETIKET <- c(SenMayo = "SenMayo (all genes)", SenMayo_immun = "SenMayo: immune-cell genes",
                 SenMayo_immun_disi = "SenMayo: stromal/adipocyte genes", Makrofaj = "Macrophage score",
                 Adiposit = "Adipocyte score", Fridman_SEN_UP = "Fridman senescence set")
zaman_etiketi <- function(acc, k) {
  k <- sub(" vs .*", "", k)
  if (acc == "GSE199063") return(c("2_year" = "2 years", "5_year" = "5 years")[k])
  if (acc == "GSE59034") return("2 years")
  if (acc == "GSE65540") return("3 months")
  if (acc == "GSE77962") return(ifelse(grepl("loss", k), "after weight loss", "after maintenance"))
  k
}
sekil_kaydet <- function(p, ad, w, h) {
  ggsave(yol(ad, ".png"), p, width = w, height = h, dpi = 300)
  tryCatch(ggsave(yol(ad, ".pdf"), p, width = w, height = h, device = cairo_pdf),
           error = function(e) ggsave(yol(ad, ".pdf"), p, width = w, height = h))
}
yayin_sekilleri <- function(sonuc, kontrol, bilesim) {
  tema <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
                                           legend.position = "bottom")
  # --- Şekil 1: Etki büyüklükleri (orman grafiği) ---
  sk <- c("SenMayo_immun", "SenMayo_immun_disi", "Makrofaj", "Adiposit")
  e <- sonuc[sonuc$grup == "tum_orneklem" & sonuc$skor %in% sk, ]
  if (nrow(e)) {
    n <- e$n_cift; d <- e$cohen_dz
    e$se <- ifelse(e$test == "eslestirilmis", sqrt(1 / n + d^2 / (2 * n)), sqrt(2 / n + d^2 / (4 * n)))
    e$kohort <- KOHORT_ETIKET[e$veri_seti]
    e$satir <- paste0(e$kohort, ", ", mapply(zaman_etiketi, e$veri_seti, e$karsilastirma))
    sira <- c(GSE65540 = 1, GSE59034 = 2, GSE199063 = 3, GSE77962 = 4)[e$veri_seti] * 10 + grepl("5_year|stable", e$karsilastirma)
    e$satir <- factor(e$satir, levels = rev(unique(e$satir[order(sira)])))
    e$skor_ad <- factor(SKOR_ETIKET[e$skor], levels = SKOR_ETIKET[sk])
    e$tur <- ifelse(e$veri_seti == "GSE77962", "Diet", "Bariatric surgery")
    p1 <- ggplot(e, aes(cohen_dz, satir, colour = tur)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      geom_errorbarh(aes(xmin = cohen_dz - 1.96 * se, xmax = cohen_dz + 1.96 * se), height = 0.25) +
      geom_point(size = 2.2) + facet_wrap(~ skor_ad, nrow = 1) +
      scale_colour_manual(values = c("Bariatric surgery" = "#2C6E9B", "Diet" = "#C8742B"), name = NULL) +
      labs(x = "Change from baseline (Cohen's d, 95% CI)", y = NULL) + tema
    sekil_kaydet(p1, "07_sekil1_etki_buyuklukleri", 13, 4.2)
  }
  # --- Şekil 2: Kontrollere göre standartlaştırılmış seyir ---
  sat <- list()
  for (acc in intersect(c("GSE199063", "GSE59034"), names(tum))) {
    f <- yol("02_skorlar_", acc, ".csv"); if (!file.exists(f)) next
    df <- read.csv(f, check.names = FALSE); a <- ELLE_AYAR[[acc]]
    z <- as.character(df$zaman); kon <- grepl("control|never|lean", z, ignore.case = TRUE)
    if (!any(kon)) next
    for (s in c("SenMayo_immun", "SenMayo_immun_disi")) {
      st <- (df[[s]] - mean(df[[s]][kon])) / sd(df[[s]][kon])
      ag <- aggregate(st, list(zaman = z), function(v) c(ort = mean(v), se = sd(v) / sqrt(length(v)), n = length(v)))
      ag <- data.frame(zaman = ag$zaman, ag$x); ag$skor <- s; ag$kohort <- KOHORT_ETIKET[acc]
      sat[[length(sat) + 1]] <- ag
    }
  }
  if (length(sat)) {
    g <- do.call(rbind, sat)
    g$asama <- ifelse(grepl("control|never|lean", g$zaman, ignore.case = TRUE), "Lean controls",
               ifelse(grepl("baseline|before", g$zaman, ignore.case = TRUE), "Pre-surgery",
               ifelse(grepl("5", g$zaman), "5 years", "2 years")))
    g$asama <- factor(g$asama, levels = c("Pre-surgery", "2 years", "5 years", "Lean controls"))
    g$skor_ad <- factor(SKOR_ETIKET[g$skor], levels = SKOR_ETIKET[c("SenMayo_immun", "SenMayo_immun_disi")])
    p2 <- ggplot(g, aes(asama, ort, colour = kohort, group = kohort)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_line(data = g[g$asama != "Lean controls", ], position = position_dodge(0.3)) +
      geom_pointrange(aes(ymin = ort - 1.96 * se, ymax = ort + 1.96 * se), position = position_dodge(0.3)) +
      facet_wrap(~ skor_ad) + scale_colour_manual(values = c("#2C6E9B", "#7A9E3B"), name = NULL) +
      labs(x = NULL, y = "Score (SD units relative to lean controls)") + tema
    sekil_kaydet(p2, "07_sekil2_kontrole_gore_seyir", 9, 4.5)
  }
  # --- Şekil 3: Doku bileşimine göre düzeltme ---
  if (!is.null(bilesim) && "ham_se" %in% colnames(bilesim)) {
    b <- bilesim[bilesim$veri_seti != "GSE77962" & bilesim$skor %in% c("SenMayo_immun", "SenMayo_immun_disi", "Fridman_SEN_UP"), ]
    if (nrow(b)) {
      u <- rbind(data.frame(b[, c("veri_seti", "karsilastirma", "skor")], tur = "Unadjusted", est = b$ham_degisim, se = b$ham_se),
                 data.frame(b[, c("veri_seti", "karsilastirma", "skor")], tur = "Adjusted for tissue composition", est = b$duzeltilmis_degisim, se = b$duzeltilmis_se))
      u$satir <- paste0(KOHORT_ETIKET[u$veri_seti], ", ", mapply(zaman_etiketi, u$veri_seti, u$karsilastirma))
      u$skor_ad <- factor(SKOR_ETIKET[u$skor], levels = SKOR_ETIKET[c("SenMayo_immun", "SenMayo_immun_disi", "Fridman_SEN_UP")])
      u$tur <- factor(u$tur, levels = c("Unadjusted", "Adjusted for tissue composition"))
      p3 <- ggplot(u, aes(est, satir, colour = tur)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
        geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), position = position_dodge(0.5)) +
        facet_wrap(~ skor_ad, scales = "free_x") + scale_colour_manual(values = c("#2C6E9B", "#B04A5A"), name = NULL) +
        labs(x = "Change in score (ssGSEA units, 95% CI)", y = NULL) + tema
      sekil_kaydet(p3, "07_sekil3_bilesim_duzeltmesi", 12, 3.8)
    }
  }
  # --- Şekil 4: Temel senesans genleri (log2 kat değişim ısı haritası) ---
  tg <- do.call(rbind, lapply(list.files(CIKTI_KLASORU, "^04_temel_senesans_genleri_", full.names = TRUE), function(f) {
    x <- read.csv(f); x$veri_seti <- sub("04_temel_senesans_genleri_(.*)\\.csv", "\\1", basename(f)); x }))
  if (!is.null(tg) && nrow(tg)) {
    tg <- tg[!grepl("never|control|lean", sub(" vs .*", "", tg$karsilastirma), ignore.case = TRUE), ]
    tg$sutun <- paste0(KOHORT_ETIKET[tg$veri_seti], "\n", mapply(zaman_etiketi, tg$veri_seti, tg$karsilastirma))
    sira <- c(GSE65540 = 1, GSE59034 = 2, GSE199063 = 3, GSE77962 = 4)[tg$veri_seti] * 10 + grepl("5_year|stable", tg$karsilastirma)
    tg$sutun <- factor(tg$sutun, levels = unique(tg$sutun[order(sira)]))
    tg$yildiz <- ifelse(tg$adj.P.Val < 0.001, "***", ifelse(tg$adj.P.Val < 0.01, "**", ifelse(tg$adj.P.Val < 0.05, "*", "")))
    lim <- max(abs(tg$logFC), na.rm = TRUE)
    p4 <- ggplot(tg, aes(sutun, gen, fill = logFC)) + geom_tile(colour = "white") +
      geom_text(aes(label = yildiz), size = 3.5) +
      scale_fill_gradient2(low = "#2C6E9B", mid = "white", high = "#B04A5A", limits = c(-lim, lim), name = "log2 FC") +
      labs(x = NULL, y = NULL, caption = "* FDR<0.05, ** FDR<0.01, *** FDR<0.001") + tema +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8), legend.position = "right")
    sekil_kaydet(p4, "07_sekil4_temel_genler", 9, 4.5)
  }
}
tryCatch(if (length(tum)) yayin_sekilleri(sonuc, kontrol, bilesim),
         error = function(e) message("!!! Şekil üretiminde hata: ", conditionMessage(e)))
# GSE65540 metadata'sını çıktıya yaz (eşleştirme kontrolü için)
for (f in list.files(CIKTI_KLASORU, "^01_metadata_GSE65540", full.names = TRUE)) cat("\n", readLines(f), sep = "\n")

message("\nBitti. Çıktılar: ", normalizePath(CIKTI_KLASORU))
message("Lütfen '", CIKTI_KLASORU, "' klasörünün tamamını (veya en azından 00_OZET.txt ve 01_metadata_*.txt) paylaşın.")
