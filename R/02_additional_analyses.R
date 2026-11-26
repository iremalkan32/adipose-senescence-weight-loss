###############################################################################
# EK ANALİZLER (hakem önerileri) — ana betik (hücre 2) çalıştıktan sonra çalıştırın.
# Ana betiğin fonksiyonlarını ve 02_skorlar_*.csv dosyalarını kullanır.
# Çıktılar: ek_analiz_cikti/ klasörü ve ekrana yazılan özet.
###############################################################################
suppressPackageStartupMessages(library(nlme))
ANA_CIKTI <- "senesans_pilot_cikti"
EK <- "ek_analiz_cikti"; dir.create(EK, showWarnings = FALSE)
ekyol <- function(...) file.path(EK, paste0(...))
SKORLAR <- c("SenMayo", "SenMayo_immun", "SenMayo_immun_disi", "Fridman_SEN_UP", "Makrofaj", "Adiposit")
kontrol_mu <- function(z) grepl("control|never|lean", z, ignore.case = TRUE)
dz_ci <- function(fark) {                       # eşleştirilmiş Cohen dz ve %95 GA (yaklaşık SE)
  n <- sum(!is.na(fark)); d <- mean(fark, na.rm = TRUE) / sd(fark, na.rm = TRUE)
  se <- sqrt(1 / n + d^2 / (2 * n)); c(n = n, dz = d, alt = d - 1.96 * se, ust = d + 1.96 * se,
                                       p = tryCatch(t.test(fark)$p.value, error = function(e) NA))
}
baslik <- function(t) cat("\n\n==================== ", t, " ====================\n", sep = "")
options(width = 200)

## ---------------------------------------------------------------------------
## A) GSE199063: karışık etkili model, 2. yıl -> 5. yıl, tam vakalar, takip kaybı
## ---------------------------------------------------------------------------
baslik("A) GSE199063 longitudinal analiz")
s <- read.csv(file.path(ANA_CIKTI, "02_skorlar_GSE199063.csv"), check.names = FALSE)
s <- s[!kontrol_mu(s$zaman), ]
s$zaman <- factor(s$zaman, levels = c("Baseline", "2_year", "5_year"))
A <- list()
for (sk in SKORLAR) {
  sd0 <- sd(s[[sk]][s$zaman == "Baseline"])
  for (ref in c("Baseline", "2_year")) {
    s$z2 <- relevel(s$zaman, ref = ref)
    m <- lme(as.formula(paste0("`", sk, "` ~ z2")), random = ~ 1 | kisi, data = s, method = "REML", na.action = na.omit)
    tt <- summary(m)$tTable
    for (r in grep("^z2", rownames(tt), value = TRUE)) {
      hedef <- sub("^z2", "", r); if (ref == "2_year" && hedef != "5_year") next
      A[[length(A) + 1]] <- data.frame(skor = sk, karsilastirma = paste(hedef, "vs", ref), model = "karisik_etkili_tum_veri",
        n_kisi = length(unique(s$kisi)), tahmin = tt[r, "Value"], std_tahmin = tt[r, "Value"] / sd0,
        alt95 = (tt[r, "Value"] - 1.96 * tt[r, "Std.Error"]) / sd0, ust95 = (tt[r, "Value"] + 1.96 * tt[r, "Std.Error"]) / sd0,
        p = tt[r, "p-value"])
    }
  }
  # Tam vakalar (üç zamanda da ölçümü olanlar)
  tam <- names(which(table(s$kisi[!is.na(s[[sk]])]) == 3))
  st <- s[s$kisi %in% tam, ]
  w <- reshape(st[, c("kisi", "zaman", sk)], idvar = "kisi", timevar = "zaman", direction = "wide")
  for (k in list(c("2_year", "Baseline"), c("5_year", "Baseline"), c("5_year", "2_year"))) {
    f <- w[[paste0(sk, ".", k[1])]] - w[[paste0(sk, ".", k[2])]]; r <- dz_ci(f)
    A[[length(A) + 1]] <- data.frame(skor = sk, karsilastirma = paste(k[1], "vs", k[2]), model = "tam_vaka_eslestirilmis_dz",
      n_kisi = r["n"], tahmin = mean(f), std_tahmin = r["dz"], alt95 = r["alt"], ust95 = r["ust"], p = r["p"])
  }
  # Takip kaybı: 5. yılda örneği olmayanların başlangıç skoru farklı mı?
  b <- s[s$zaman == "Baseline", ]; var5 <- b$kisi %in% s$kisi[s$zaman == "5_year"]
  if (sum(!var5) >= 3) {
    tt <- t.test(b[[sk]][!var5], b[[sk]][var5])
    dd <- (mean(b[[sk]][!var5]) - mean(b[[sk]][var5])) / sqrt((var(b[[sk]][!var5]) + var(b[[sk]][var5])) / 2)
    A[[length(A) + 1]] <- data.frame(skor = sk, karsilastirma = paste0("baslangic: 5y eksik (n=", sum(!var5), ") vs 5y var (n=", sum(var5), ")"),
      model = "takip_kaybi_Welch", n_kisi = nrow(b), tahmin = diff(rev(tt$estimate)), std_tahmin = dd,
      alt95 = NA, ust95 = NA, p = tt$p.value)
  }
}
A <- do.call(rbind, A); rownames(A) <- NULL
write.csv(A, ekyol("A_GSE199063_longitudinal.csv"), row.names = FALSE); print(A, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------------------
## B) Kontrollerle eşdeğerlik (TOST), sınırlar ±0.2 ve ±0.5 SD
## ---------------------------------------------------------------------------
baslik("B) Zayıf kontrollerle eşdeğerlik (TOST)")
B <- list()
for (acc in c("GSE199063", "GSE59034")) {
  s <- read.csv(file.path(ANA_CIKTI, paste0("02_skorlar_", acc, ".csv")), check.names = FALSE)
  kon <- kontrol_mu(s$zaman); zl <- setdiff(unique(s$zaman), unique(s$zaman[kon]))
  for (z in zl) for (sk in SKORLAR) {
    x <- s[[sk]][s$zaman == z]; y <- s[[sk]][kon]
    sp <- sqrt((var(x) + var(y)) / 2); d <- (mean(x) - mean(y)) / sp
    se <- sqrt(1 / length(x) + 1 / length(y) + d^2 / (2 * (length(x) + length(y))))
    satir <- data.frame(veri_seti = acc, zaman = z, skor = sk, n = length(x), n_kontrol = length(y), d = d,
                        ga90_alt = d - 1.645 * se, ga90_ust = d + 1.645 * se, p_fark = t.test(x, y)$p.value)
    for (sinir in c(0.2, 0.5)) {
      dl <- sinir * sp
      p1 <- t.test(x, y, mu = -dl, alternative = "greater")$p.value
      p2 <- t.test(x, y, mu = dl, alternative = "less")$p.value
      satir[[paste0("p_TOST_", sinir)]] <- max(p1, p2)
    }
    B[[length(B) + 1]] <- satir
  }
}
B <- do.call(rbind, B); write.csv(B, ekyol("B_TOST_kontrol.csv"), row.names = FALSE); print(B, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------------------
## Verileri yeniden yükle (indirilen dosyalar önbellekte)
## ---------------------------------------------------------------------------
eslestir <- function(d) {     # her takip zamanı için başlangıçla eşleşen örnek indeksleri
  a <- sutunlari_algila(d$acc, d$pd); pd <- if (!is.null(a$pd)) a$pd else d$pd
  bul <- function(ad) { if (is.null(ad) || ad %in% colnames(pd)) return(ad); b <- sub(":ch1$", "", ad); if (b %in% colnames(pd)) b else make.names(ad) }
  z <- as.character(pd[[bul(a$zaman)]]); k <- as.character(pd[[bul(a$kisi)]])
  tak <- setdiff(unique(z), c(a$baslangic, a$kontrol))
  out <- lapply(tak, function(t) { o <- intersect(k[z == a$baslangic], k[z == t])
    list(zaman = t, s0 = which(z == a$baslangic)[match(o, k[z == a$baslangic])], s1 = which(z == t)[match(o, k[z == t])]) })
  names(out) <- tak; list(es = out, zaman = z, kisi = k, kontrol = a$kontrol, baslangic = a$baslangic)
}
D <- list()
for (acc in c("GSE199063", "GSE59034", "GSE77962", "GSE65540"))
  D[[acc]] <- tryCatch(if (acc %in% RNASEQ_SETLERI) rnaseq_yukle(acc) else gse_yukle(acc),
                       error = function(e) { message("!!! ", acc, ": ", conditionMessage(e)); NULL })

## ---------------------------------------------------------------------------
## C) İfade düzeyine göre eşleştirilmiş rastgele gen setleri, SASP dışı SenMayo,
##    SenMayo ile örtüşmeyen Fridman seti, SenMayo–Fridman korelasyonu
## ---------------------------------------------------------------------------
baslik("C) Eşleştirilmiş rastgele setler ve alt setler")
SASP_DESENI <- "^(CCL|CXCL|IL[0-9]|CSF[0-9]|TNF|MMP|SERPIN|IGFBP|WNT|BMP|FGF[0-9]|VEGF|PGF$|HGF$|ANG$|ANGPT|AREG|EREG|GDF|INHA|LIF$|TIMP|PLAU|SPP1|CTSB|KITLG|NRG|SEMA|HMGB|ICAM|EDN|ESM1|PAPPA|IGF[0-9]|CXCL|C3$|CD55|EGF$|FAS$|GEM$|JUN$|MIF$|PTGES|PTBP1|TGFB)"
sasp_disi <- setler$SenMayo[!grepl(SASP_DESENI, setler$SenMayo)]
fr_ozel <- setdiff(setler$Fridman_SEN_UP, setler$SenMayo)
cat("SenMayo–Fridman ortak gen sayısı:", length(intersect(setler$SenMayo, setler$Fridman_SEN_UP)),
    "| SASP dışı SenMayo gen sayısı:", length(sasp_disi), "| SenMayo ile örtüşmeyen Fridman:", length(fr_ozel), "\n")
writeLines(sasp_disi, ekyol("C_SenMayo_SASP_disi_genler.txt"))
set.seed(2026); C <- list()
for (acc in names(D)) {
  d <- D[[acc]]; if (is.null(d)) next
  e <- eslestir(d); g_sen <- intersect(setler$SenMayo, rownames(d$x))
  ort <- rowMeans(d$x); kutu <- cut(rank(ort), 20, labels = FALSE); names(kutu) <- rownames(d$x)
  aday <- setdiff(rownames(d$x), setler$SenMayo)
  havuz <- split(aday, kutu[aday])
  rastgele_set <- function() unlist(lapply(g_sen, function(g) sample(havuz[[as.character(kutu[g])]], 1)))
  sk <- ssgsea_skor(d$x, list(SenMayo = g_sen, SASP_disi = sasp_disi, Fridman_ozel = fr_ozel, Fridman = setler$Fridman_SEN_UP))
  for (t in names(e$es)) {
    s0 <- e$es[[t]]$s0; s1 <- e$es[[t]]$s1; if (length(s0) < 3) next
    degisim <- function(genler) { z <- z_skor(d$x, genler); mean(z[s1] - z[s0]) }
    gercek <- degisim(g_sen); rast <- replicate(1000, degisim(rastgele_set()))
    satir <- data.frame(veri_seti = acc, karsilastirma = paste(t, "vs", e$baslangic), n = length(s0),
      gercek_degisim = gercek, eslesmis_rastgele_ort = mean(rast),
      p_eslesmis_rastgele = (sum(abs(rast - mean(rast)) >= abs(gercek - mean(rast))) + 1) / 1001)
    for (alt in c("SASP_disi", "Fridman_ozel")) { r <- dz_ci(sk[alt, s1] - sk[alt, s0])
      satir[[paste0(alt, "_dz")]] <- r["dz"]; satir[[paste0(alt, "_p")]] <- r["p"] }
    satir$SenMayo_Fridman_r <- cor(sk["SenMayo", ], sk["Fridman", ])
    C[[length(C) + 1]] <- satir
  }
}
C <- do.call(rbind, C); write.csv(C, ekyol("C_eslesmis_rastgele_ve_alt_setler.csv"), row.names = FALSE); print(C, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------------------
## D) GSE59034 ile GSE199063 arasında katılımcı örtüşmesi: genotipe bağlı
##    ("kişiye özgü") ifade parmak izi. Her iki veri setinde de kişi içi tekrar
##    ölçümlerde çok kararlı (yüksek ICC) genler seçilir; kişiler arası korelasyon
##    aynı kişi (pozitif kontrol) ve farklı kişi (boş dağılım) ile karşılaştırılır.
## ---------------------------------------------------------------------------
baslik("D) Katılımcı örtüşmesi (ifade parmak izi)")
icc <- function(x, kisi) {    # tek yönlü ANOVA ICC(1), yalnızca >=2 örneği olan kişiler
  tut <- kisi %in% names(which(table(kisi) >= 2)); x <- x[, tut, drop = FALSE]; k <- factor(kisi[tut])
  n0 <- mean(table(k)); m <- t(apply(x, 1, function(v) tapply(v, k, mean)))
  msb <- apply(m, 1, var) * n0; msw <- rowMeans((x - m[, as.integer(k)])^2) * ncol(x) / (ncol(x) - nlevels(k))
  (msb - msw) / (msb + (n0 - 1) * msw)
}
if (!is.null(D$GSE199063) && !is.null(D$GSE59034)) {
  a <- D$GSE199063; b <- D$GSE59034; ea <- eslestir(a); eb <- eslestir(b)
  ortak <- intersect(rownames(a$x), rownames(b$x))
  ia <- icc(a$x[ortak, ], ea$kisi); ib <- icc(b$x[ortak, ], eb$kisi)
  sdmin <- pmin(apply(a$x[ortak, ], 1, sd), apply(b$x[ortak, ], 1, sd))
  aday <- ortak[ia > 0.8 & ib > 0.8 & sdmin > 0.3]
  aday <- aday[order(-pmin(ia[aday], ib[aday]))][1:min(300, length(aday))]
  cat("Parmak izi gen sayısı:", length(aday), "| ilk 20:", paste(head(aday, 20), collapse = ", "), "\n")
  za <- t(scale(t(a$x[aday, ]))); zb <- t(scale(t(b$x[aday, ])))
  pa <- sapply(split(seq_along(ea$kisi), ea$kisi), function(i) rowMeans(za[, i, drop = FALSE]))
  pb <- sapply(split(seq_along(eb$kisi), eb$kisi), function(i) rowMeans(zb[, i, drop = FALSE]))
  # Pozitif kontrol: GSE199063'te aynı kişinin başlangıç ve 2. yıl örnekleri (örnek düzeyinde)
  s0 <- ea$es[["2_year"]]$s0; s1 <- ea$es[["2_year"]]$s1
  ayni <- sapply(seq_along(s0), function(i) cor(za[, s0[i]], za[, s1[i]]))
  farkli <- as.vector(cor(za[, s0], za[, s1])[upper.tri(matrix(0, length(s0), length(s0)))])
  cat(sprintf("GSE199063 aynı kişi r: medyan %.2f (min %.2f) | farklı kişi r: medyan %.2f, %%99.9 = %.2f, maks %.2f\n",
              median(ayni), min(ayni), median(farkli), quantile(farkli, 0.999), max(farkli)))
  capraz <- cor(pb, pa)        # satır: GSE59034 kişileri, sütun: GSE199063 kişileri
  esik <- quantile(farkli, 0.999)
  en_iyi <- data.frame(GSE59034_kisi = rownames(capraz), en_iyi_GSE199063_kisi = colnames(capraz)[apply(capraz, 1, which.max)],
                       r_en_iyi = apply(capraz, 1, max), r_ikinci = apply(capraz, 1, function(v) sort(v, decreasing = TRUE)[2]))
  en_iyi$karsilikli_en_iyi <- sapply(seq_len(nrow(en_iyi)), function(i)
    rownames(capraz)[which.max(capraz[, en_iyi$en_iyi_GSE199063_kisi[i]])] == en_iyi$GSE59034_kisi[i])
  esik <- max(farkli)          # farklı kişiler arasında görülen en yüksek korelasyon
  cat(sprintf("Eşik: r > %.2f (farklı kişi maksimumu), karşılıklı en iyi eşleşme ve ikinci adaydan en az 0.15 yüksek\n", esik))
  en_iyi$olasi_ayni_kisi <- en_iyi$r_en_iyi > esik & en_iyi$karsilikli_en_iyi & (en_iyi$r_en_iyi - en_iyi$r_ikinci) > 0.15
  print(en_iyi[order(-en_iyi$r_en_iyi), ], digits = 3, row.names = FALSE)
  cat("\nOlası aynı kişi sayısı (yukarıdaki üç ölçüt):", sum(en_iyi$olasi_ayni_kisi),
      "/", nrow(en_iyi), "\n")
  write.csv(en_iyi, ekyol("D_katilimci_ortusmesi.csv"), row.names = FALSE)
  cat("\nGSE59034 örnek bilgileri (ilk 40 satır):\n")
  print(head(b$pd[, c("title", grep(":ch1$", colnames(b$pd), value = TRUE))], 40))
  cat("\nGSE59034 seri bilgileri:\n")
  f <- list.files(INDIRME_KLASORU, "GSE59034.*series_matrix", full.names = TRUE)[1]
  if (!is.na(f)) cat(grep("^!Series_(title|summary|overall_design|contributor|pubmed_id|submission_date|contact_name|contact_institute|contact_country)",
                          readLines(gzfile(f), n = 300), value = TRUE), sep = "\n")
}

## ---------------------------------------------------------------------------
## E) Duyarlılık: RNA-seq için TMM; mikroarray için prob ortalaması
## ---------------------------------------------------------------------------
baslik("E) Normalizasyon ve prob birleştirme duyarlılığı")
yedek <- list(CIKTI_KLASORU, DUYARLILIK, RASTGELE_SET_SAYISI, PROB_BIRLESTIRME)
DUYARLILIK <- FALSE; RASTGELE_SET_SAYISI <- 10
E <- list()
calistir <- function(acc, etiket) {
  CIKTI_KLASORU <<- file.path(EK, etiket); dir.create(CIKTI_KLASORU, showWarnings = FALSE)
  r <- tryCatch(analiz_et(if (acc %in% RNASEQ_SETLERI) rnaseq_yukle(acc) else gse_yukle(acc), setler),
                error = function(e) { message("!!! ", acc, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  x <- r$sonuc[r$sonuc$grup == "tum_orneklem" & r$sonuc$skor %in% SKORLAR, c("veri_seti", "karsilastirma", "skor", "n_cift", "cohen_dz", "p_eslestirilmis_t")]
  x$duyarlilik <- etiket; x
}
RNASEQ_TMM <- TRUE;  E[[1]] <- calistir("GSE65540", "TMM"); RNASEQ_TMM <- FALSE
PROB_BIRLESTIRME <- "ortalama"
for (acc in c("GSE199063", "GSE59034", "GSE77962")) E[[length(E) + 1]] <- calistir(acc, "prob_ortalamasi")
CIKTI_KLASORU <- yedek[[1]]; DUYARLILIK <- yedek[[2]]; RASTGELE_SET_SAYISI <- yedek[[3]]; PROB_BIRLESTIRME <- yedek[[4]]
E <- do.call(rbind, E); write.csv(E, ekyol("E_normalizasyon_prob_duyarlilik.csv"), row.names = FALSE); print(E, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------------------
## F) İkincil testler için Benjamini–Hochberg FDR
## ---------------------------------------------------------------------------
baslik("F) İkincil testler için BH-FDR")
S <- read.csv(file.path(ANA_CIKTI, "03_eslestirilmis_testler.csv"))
S <- S[S$grup == "tum_orneklem", ]
ikincil <- S$skor != "SenMayo"
S$FDR_ikincil <- NA; S$FDR_ikincil[ikincil] <- p.adjust(S$p_eslestirilmis_t[ikincil], "BH")
birincil <- S$skor == "SenMayo"
S$Holm_birincil <- NA; S$Holm_birincil[birincil] <- p.adjust(S$p_eslestirilmis_t[birincil], "holm")
S2 <- S[S$skor %in% SKORLAR, c("veri_seti", "karsilastirma", "skor", "n_cift", "cohen_dz", "p_eslestirilmis_t", "FDR_ikincil", "Holm_birincil")]
write.csv(S2, ekyol("F_FDR.csv"), row.names = FALSE); print(S2, digits = 3, row.names = FALSE)
K <- tryCatch(read.csv(file.path(ANA_CIKTI, "03_bilesim_duzeltmesi.csv")), error = function(e) NULL)
if (!is.null(K)) { K$FDR <- p.adjust(K$duzeltilmis_p, "BH")
  cat("\nBileşim düzeltmesi, BH-FDR:\n"); print(K[, c("veri_seti", "karsilastirma", "skor", "duzeltilmis_degisim", "duzeltilmis_p", "FDR")], digits = 3, row.names = FALSE)
  write.csv(K, ekyol("F_bilesim_FDR.csv"), row.names = FALSE) }

cat("\n\nEK ANALİZLER BİTTİ. Klasör: ", normalizePath(EK), "\n")
