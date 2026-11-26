###############################################################################
# D2) Katılımcı eşleştirmesi: kişiye özgü ifade parmak izi
#  1) GSE199063'te (gerçek kişi kodları var) kişi içinde çok kararlı genler seçilir.
#  2) Yöntem GSE199063 içinde doğrulanır (aynı kişi vs farklı kişi).
#  3) GSE59034 içinde "rep" numarasına göre yapılan önce/sonra eşleştirmesi kontrol edilir.
#  4) GSE59034 örnekleri GSE199063 kişileriyle karşılaştırılır (ortak katılımcı var mı?).
###############################################################################
options(width = 220)
dir.create("ek_analiz_cikti", showWarnings = FALSE)
a <- gse_yukle("GSE199063"); b <- gse_yukle("GSE59034")
pa <- a$pd; pb <- b$pd
kisiA <- as.character(pa[["subject_id:ch1"]]); zamA <- as.character(pa[["time_point:ch1"]])
zamB <- sub("obesity status: ", "", as.character(pb[["obesity status:ch1"]]))
if (all(is.na(zamB)) || !length(zamB)) zamB <- ifelse(grepl("before", pb$title), "before", ifelse(grepl("after", pb$title), "after", "never"))
zamB <- ifelse(grepl("before", zamB), "before", ifelse(grepl("after", zamB), "after", "never"))
repB <- sub(".*rep", "", pb$title)
cat("GSE59034 grupları:", paste(names(table(zamB)), table(zamB), collapse = ", "), "\n")

ortak <- intersect(rownames(a$x), rownames(b$x))
xa <- a$x[ortak, ]; xb <- b$x[ortak, ]
tamam <- rowSums(is.na(xa)) == 0 & rowSums(is.na(xb)) == 0
xa <- xa[tamam, ]; xb <- xb[tamam, ]
cat("Ortak gen sayısı:", nrow(xa), "\n")

# 1) GSE199063 kişi içi kararlılık (ICC), yalnızca ameliyatlı kadınlar (>=2 örnek)
icc <- function(x, kisi) {
  tut <- kisi %in% names(which(table(kisi) >= 2)); x <- x[, tut, drop = FALSE]; k <- factor(kisi[tut])
  n0 <- mean(table(k)); m <- t(apply(x, 1, function(v) tapply(v, k, mean)))
  msb <- apply(m, 1, var) * n0; msw <- rowSums((x - m[, as.integer(k)])^2) / (ncol(x) - nlevels(k))
  (msb - msw) / (msb + (n0 - 1) * msw)
}
# Zaman etkisini çıkar: her gen için zaman noktası ortalaması çıkarılır (ameliyatın genel etkisi kişi içi varyansı şişirmesin)
merkez <- function(x, grup) { for (g in unique(grup)) { i <- grup == g; x[, i] <- x[, i] - rowMeans(x[, i, drop = FALSE]) }; x }
xa <- merkez(xa, zamA)
xb <- merkez(xb, zamB)
ia <- icc(xa, kisiA)
cat("GSE199063 ICC dağılımı:\n"); print(summary(ia)); print(table(cut(ia, c(-Inf, .5, .7, .8, .9, Inf))))
sdA <- apply(xa, 1, sd); sdB <- apply(xb, 1, sd)
aday <- names(ia)[!is.na(ia) & ia > 0.5 & sdA > 0.15 & sdB > 0.15]
aday <- aday[order(-ia[aday])][1:min(200, length(aday))]
cat("Seçilen genlerin ICC aralığı:", round(range(ia[aday]), 2), "\n")
cat("Parmak izi gen sayısı:", length(aday), "| ilk 30:", paste(head(aday, 30), collapse = ", "), "\n")

za <- t(scale(t(xa[aday, ]))); zb <- t(scale(t(xb[aday, ])))

# 2) Doğrulama GSE199063 içinde
b0 <- which(zamA == "Baseline"); y2 <- which(zamA == "2_year")
o <- intersect(kisiA[b0], kisiA[y2]); i0 <- b0[match(o, kisiA[b0])]; i2 <- y2[match(o, kisiA[y2])]
R <- cor(za[, i0], za[, i2])
ayni <- diag(R); farkli <- R[row(R) != col(R)]
cat(sprintf("\nGSE199063 doğrulama: aynı kişi r medyan %.2f (min %.2f) | farklı kişi r medyan %.2f, maks %.2f | doğru eşleşme (en iyi aday aynı kişi): %d/%d\n",
            median(ayni), min(ayni), median(farkli), max(farkli), sum(apply(R, 1, which.max) == seq_len(nrow(R))), nrow(R)))
esik <- as.numeric(quantile(farkli, 0.999))
cat(sprintf("Farklı kişi %%99.9 = %.2f; aynı kişilerin %%%d'i bu eşiğin üstünde\n", esik, round(100 * mean(ayni > esik))))

# 3) GSE59034 içinde önce/sonra eşleştirmesi
ib <- which(zamB == "before"); iaf <- which(zamB == "after"); inv <- which(zamB == "never")
RB <- cor(zb[, ib], zb[, iaf]); dimnames(RB) <- list(paste0("before_rep", repB[ib]), paste0("after_rep", repB[iaf]))
en_iyi <- apply(RB, 1, which.max)
tab3 <- data.frame(once = rownames(RB), en_iyi_sonra = colnames(RB)[en_iyi], r_en_iyi = apply(RB, 1, max),
                   r_ikinci = apply(RB, 1, function(v) sort(v, decreasing = TRUE)[2]),
                   rep_eslesmesi_r = sapply(seq_len(nrow(RB)), function(i) { j <- match(repB[ib][i], repB[iaf]); if (is.na(j)) NA else RB[i, j] }))
tab3$rep_numarasi_ayni_kisi_mi <- sub("before_", "", tab3$once) == sub("after_", "", tab3$en_iyi_sonra)
cat("\nGSE59034: her 'before' örneğinin en çok benzediği 'after' örneği (rep numarası eşleştirmesi doğru mu?)\n")
print(tab3, digits = 3, row.names = FALSE)
cat("rep numarası doğru eşleşme sayısı:", sum(tab3$rep_numarasi_ayni_kisi_mi), "/", nrow(tab3), "\n")
write.csv(tab3, "ek_analiz_cikti/D3_GSE59034_ic_eslesme.csv", row.names = FALSE)
RN <- cor(zb[, inv], zb[, c(ib, iaf)])
cat(sprintf("Kontrollerin ameliyatlı örneklere en yüksek benzerliği: maks r = %.2f\n", max(RN)))

# 4) GSE59034 örnekleri vs GSE199063 kişileri (kişi ortalaması parmak izi)
kA <- split(seq_along(kisiA), kisiA)
profA <- sapply(kA, function(i) rowMeans(za[, i, drop = FALSE]))
zamA_kisi <- sapply(kA, function(i) if (any(zamA[i] == "Control")) "Control" else "RYGB")
C <- cor(zb, profA)                         # satır: GSE59034 örnekleri, sütun: GSE199063 kişileri
sonuc <- data.frame(GSE59034_ornek = paste0(zamB, "_rep", repB), en_iyi_GSE199063 = colnames(C)[apply(C, 1, which.max)],
                    grup_GSE199063 = zamA_kisi[colnames(C)[apply(C, 1, which.max)]],
                    r_en_iyi = apply(C, 1, max), r_ikinci = apply(C, 1, function(v) sort(v, decreasing = TRUE)[2]))
sonuc$fark <- sonuc$r_en_iyi - sonuc$r_ikinci
sonuc$olasi_ayni_kisi <- sonuc$r_en_iyi > esik & sonuc$fark > 0.15
cat(sprintf("\nGSE59034 örnekleri ile GSE199063 kişileri (eşik: r > %.2f ve ikinci adaydan en az 0.15 fazla):\n", esik))
print(sonuc[order(-sonuc$r_en_iyi), ], digits = 3, row.names = FALSE)
cat("\nOlası ortak katılımcı içeren GSE59034 örnek sayısı:", sum(sonuc$olasi_ayni_kisi), "/", nrow(sonuc), "\n")
cat("Eşleşen farklı GSE199063 kişisi sayısı:", length(unique(sonuc$en_iyi_GSE199063[sonuc$olasi_ayni_kisi])), "\n")
write.csv(sonuc, "ek_analiz_cikti/D3_capraz_eslesme.csv", row.names = FALSE)
cat("\nD3 bitti.\n")
