suppressPackageStartupMessages(library(ggplot2))
cikti <- "figures"
dir.create(cikti, showWarnings = FALSE, recursive = TRUE)

kaydet <- function(p, ad, w, h) {
  ggsave(file.path(cikti, paste0(ad, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(cikti, paste0(ad, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(cikti, paste0(ad, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}
tema <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey92"), legend.position = "bottom",
        text = element_text(colour = "black"), axis.text = element_text(colour = "black"))

skor_ad <- c(imm = "SenMayo: immune-cell genes", str = "SenMayo: stromal/adipocyte genes",
             mac = "Macrophage score", adi = "Adipocyte score", fri = "Fridman senescence set")

# ---------------- Figure 1: paired effect sizes ----------------
f1 <- read.table(header = TRUE, text = "
satir n imm str mac adi tur
'Surgery (GSE65540), 3 months' 22 -0.03175 0.06068 0.132 1.2516 S
'Surgery (GSE59034), 2 years' 16 -1.69069 -2.12729 -2.38362 0.83056 S
'Surgery (GSE199063), 2 years' 49 -1.35406 -1.38248 -1.91062 1.31454 S
'Surgery (GSE199063), 5 years' 38 -1.70144 -1.78132 -1.76362 1.88203 S
'Diet (GSE77962), after weight loss' 48 0.26646 0.23153 0.4147 -0.13615 D
'Diet (GSE77962), after maintenance' 49 0.05177 0.02435 -0.05433 -0.00847 D")
e <- do.call(rbind, lapply(c("imm", "str", "mac", "adi"), function(s)
  data.frame(satir = f1$satir, n = f1$n, d = f1[[s]], skor = skor_ad[s], tur = f1$tur)))
e$se <- sqrt(1 / e$n + e$d^2 / (2 * e$n))
e$satir <- factor(e$satir, levels = rev(f1$satir))
e$skor <- factor(e$skor, levels = skor_ad[c("imm", "str", "mac", "adi")])
e$tur <- factor(ifelse(e$tur == "S", "Bariatric surgery", "Diet"), levels = c("Bariatric surgery", "Diet"))
p1 <- ggplot(e, aes(d, satir, shape = tur)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = d - 1.96 * se, xmax = d + 1.96 * se), width = 0.25, orientation = "y") +
  geom_point(size = 2.4, fill = "white", stroke = 0.8) + facet_wrap(~ skor, nrow = 1) +
  scale_shape_manual(values = c("Bariatric surgery" = 16, "Diet" = 24), name = NULL) +
  labs(x = "Change from baseline (paired Cohen's dz, 95% CI)", y = NULL) + tema
kaydet(p1, "Figure1", 13, 4.2)

# ---------------- Figure 2: distance from lean controls ----------------
f2 <- read.table(header = TRUE, text = "
kohort asama n nk imm str
'Surgery (GSE199063)' Pre-surgery 50 28 1.2041 1.6226
'Surgery (GSE199063)' '2 years' 49 28 -0.4533 -0.222
'Surgery (GSE199063)' '5 years' 38 28 -0.8836 -0.5311
'Surgery (GSE59034)' Pre-surgery 16 16 1.4379 2.2916
'Surgery (GSE59034)' '2 years' 16 16 -0.2153 0.3832")
g <- do.call(rbind, lapply(c("imm", "str"), function(s)
  data.frame(f2[, 1:4], d = f2[[s]], skor = skor_ad[s])))
g$se <- sqrt((g$n + g$nk) / (g$n * g$nk) + g$d^2 / (2 * (g$n + g$nk)))
g$asama <- factor(g$asama, levels = c("Pre-surgery", "2 years", "5 years"))
g$skor <- factor(g$skor, levels = skor_ad[c("imm", "str")])
dg <- position_dodge(0.3)
p2 <- ggplot(g, aes(asama, d, shape = kohort, linetype = kohort, group = kohort)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(position = dg) +
  geom_errorbar(aes(ymin = d - 1.96 * se, ymax = d + 1.96 * se), width = 0.15, position = dg, linetype = 1) +
  geom_point(size = 2.6, fill = "white", stroke = 0.8, position = dg) +
  facet_wrap(~ skor) +
  scale_shape_manual(values = c(16, 22), name = NULL) + scale_linetype_manual(values = c(1, 2), name = NULL) +
  labs(x = NULL, y = "Difference from lean controls\n(Cohen's d, 95% CI; 0 = lean level)") + tema
kaydet(p2, "Figure2", 9, 4.5)

# ---------------- Figure 3: composition adjustment ----------------
f3 <- read.table(header = TRUE, text = "
satir skor ham ham_se duz duz_se
'Surgery (GSE65540), 3 months' imm -0.000747 0.00502 -0.001265 0.00464
'Surgery (GSE65540), 3 months' str 0.00092 0.00323 0.005877 0.00529
'Surgery (GSE65540), 3 months' fri 3.14e-05 0.00196 0.001769 0.00254
'Surgery (GSE59034), 2 years' imm -0.0639 0.00945 -0.004177 0.02745
'Surgery (GSE59034), 2 years' str -0.0367 0.00431 -0.025422 0.01476
'Surgery (GSE59034), 2 years' fri -0.021 0.00216 -0.021278 0.0076
'Surgery (GSE199063), 2 years' imm -0.047 0.00495 0.009606 0.01007
'Surgery (GSE199063), 2 years' str -0.0416 0.0043 -0.023399 0.01073
'Surgery (GSE199063), 2 years' fri -0.0185 0.00238 -0.010212 0.00621
'Surgery (GSE199063), 5 years' imm -0.0643 0.00613 -0.011314 0.01331
'Surgery (GSE199063), 5 years' str -0.05 0.00455 -0.017574 0.01076
'Surgery (GSE199063), 5 years' fri -0.0217 0.00241 -0.013706 0.00646")
u <- rbind(data.frame(satir = f3$satir, skor = f3$skor, tur = "Unadjusted", est = f3$ham, se = f3$ham_se),
           data.frame(satir = f3$satir, skor = f3$skor, tur = "Adjusted for tissue composition", est = f3$duz, se = f3$duz_se))
u$satir <- factor(u$satir, levels = rev(unique(f3$satir)))
u$skor <- factor(skor_ad[u$skor], levels = skor_ad[c("imm", "str", "fri")])
u$tur <- factor(u$tur, levels = c("Unadjusted", "Adjusted for tissue composition"))
dg3 <- position_dodge(0.55)
p3 <- ggplot(u, aes(est, satir, shape = tur, group = tur)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), width = 0.25, orientation = "y", position = dg3) +
  geom_point(size = 2.4, fill = "white", stroke = 0.8, position = dg3) +
  facet_wrap(~ skor, scales = "free_x") +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  labs(x = "Change in score (ssGSEA units, 95% CI)", y = NULL) + tema
kaydet(p3, "Figure3", 12, 3.8)

# ---------------- Figure 4: classical markers ----------------
genler <- c("CDKN1A", "CDKN2A", "CDKN2B", "GLB1", "LMNB1", "MKI67", "SERPINE1", "TP53")
sut <- list(
  "Surgery (GSE65540)\n3 months" = list(fc = c(0.01796, -0.00761, -0.47716, -0.33585, -0.0726, -0.2526, -0.01174, 0.08941),
                                        q = c(0.902, 0.9695, 0.00498, 0.00092, 0.66144, 0.47147, 0.97688, 0.10862)),
  "Surgery (GSE59034)\n2 years" = list(fc = c(-0.01313, 0.00562, -0.98375, -1.03188, -0.1725, -0.68625, -1.08375, -0.21875),
                                       q = c(9.52e-01, 9.63e-01, 4.16e-08, 2.62e-09, 4.24e-02, 1.28e-04, 2.08e-08, 7.94e-04)),
  "Surgery (GSE199063)\n2 years" = list(fc = c(-0.476, -0.182, -1.03, -1.821, -0.554, -0.87, -1.63, -0.366),
                                        q = c(1.15e-01, 3.98e-03, 2.70e-12, 6.27e-23, 9.95e-07, 4.18e-10, 1.28e-16, 9.75e-05)),
  "Surgery (GSE199063)\n5 years" = list(fc = c(-0.561, -0.144, -0.901, -1.464, -0.839, -0.876, -1.492, -0.188),
                                        q = c(8.82e-02, 6.31e-02, 1.12e-08, 7.33e-16, 7.49e-11, 9.63e-09, 6.84e-13, 1.39e-01)),
  "Diet (GSE77962)\nafter weight loss" = list(fc = c(0.0983, -0.0302, -0.2224, -0.0593, 0.1173, 0.1778, -0.0182, 0.0803),
                                              q = c(0.16508, 0.75978, 0.00043, 0.4673, 0.2017, 0.1931, 0.90025, 0.14544)),
  "Diet (GSE77962)\nafter maintenance" = list(fc = c(-0.0183, -0.0378, -0.2193, -0.1369, 0.0332, -0.046, -0.2399, 0.0198),
                                              q = c(0.88986, 0.79775, 0.0015, 0.12568, 0.83711, 0.85058, 0.07328, 0.84058)))
tg <- do.call(rbind, lapply(names(sut), function(k) data.frame(sutun = k, gen = genler, fc = sut[[k]]$fc, q = sut[[k]]$q)))
tg$sutun <- factor(tg$sutun, levels = names(sut))
tg$gen <- factor(tg$gen, levels = rev(genler))
tg$yildiz <- ifelse(tg$q < 0.001, "***", ifelse(tg$q < 0.01, "**", ifelse(tg$q < 0.05, "*", "")))
tg$etiket <- paste0(sprintf("%.2f", tg$fc), tg$yildiz)
tg$yazi <- ifelse(tg$fc < -1, "white", "black")
p4 <- ggplot(tg, aes(sutun, gen, fill = fc)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = etiket, colour = yazi), size = 3) + scale_colour_identity() +
  scale_fill_gradient2(low = "grey15", mid = "white", high = "grey55", midpoint = 0, limits = c(-2, 2),
                       name = "log2 fold\nchange") +
  labs(x = NULL, y = NULL, caption = "Values are log2 fold changes. * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001") +
  tema + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(face = "italic"),
               legend.position = "right", panel.grid = element_blank())
kaydet(p4, "Figure4", 10, 4.8)
cat("done\n")
