# Charts for the VaR backtest, saved as PNG files in figures/.
# Base R graphics only, so no extra packages are needed.
# Expects the forecast table with breach flags (from add_breach_flags()) and
# the per-year table (from breaches_by_year()).

source(file.path("R", "config.R"))

FIG_DIR <- "figures"

# Colours: one fixed categorical slot per model (validated as a set for
# colour-vision-deficiency separation), plus neutral greys for returns, grid
# and text. A model keeps its colour in every chart. Identity is also carried
# by marker shape, legends and count labels, never by colour alone.
# The order is the validated slot order and also the bar order in the
# per-year chart: normal-based models first, then the Student-t ones. The
# zoom chart uses the first three slots, which also pass when every pair is
# compared (not only neighbours).
MODEL_COL <- c(historical  = "#2a78d6",   # blue
               normal      = "#eb6834",   # orange
               ewma_normal = "#1baf7a",   # aqua
               student_t   = "#eda100",   # yellow
               ewma_t      = "#e87ba4")   # magenta
MODEL_LABEL <- c(historical  = "Historical",
                 normal      = "Normal",
                 student_t   = "Student-t (df 5)",
                 ewma_normal = "EWMA normal",
                 ewma_t      = "EWMA Student-t (df 5)")
COL <- list(
  hist    = MODEL_COL[["historical"]],
  normal  = MODEL_COL[["normal"]],
  ewma    = MODEL_COL[["ewma_normal"]],
  returns = "#a3a29a",
  grid    = "#e1e0d9",
  axis    = "#c3c2b7",
  text    = "#0b0b0b",
  text2   = "#52514e",
  bg      = "#fcfcfb"
)
PCH <- list(hist = 17, normal = 1, ewma = 18)   # triangle / open circle / diamond

open_png <- function(file, width = 2400, height = 1300, bottom = 6.2) {
  dir.create(FIG_DIR, showWarnings = FALSE)
  png(file.path(FIG_DIR, file), width = width, height = height, res = 200, bg = COL$bg)
  par(mar = c(bottom, 5, 4.6, 1.5), family = "sans", col.axis = COL$text2,
      col.lab = COL$text, fg = COL$axis, las = 1, mgp = c(3.4, 0.7, 0))
}

pct_axis <- function(side = 2) {
  at <- pretty(par("usr")[3:4])
  axis(side, at = at, labels = paste0(at, "%"), lwd = 0, lwd.ticks = 0)
  abline(h = at, col = COL$grid, lwd = 0.8)
}

# Horizontal legend under the plot, each entry only as wide as its own text.
bottom_legend <- function(labels, cex, inset = -0.2, x.intersp = 0.5, ...) {
  par(xpd = NA)
  legend("bottom", inset = c(0, inset), horiz = TRUE, bty = "n", cex = cex,
         x.intersp = x.intersp,
         text.col = COL$text, legend = labels,
         text.width = strwidth(labels, cex = cex) + strwidth("MMM", cex = cex), ...)
}

title_block <- function(main, sub) {
  mtext(main, side = 3, line = 2.6, adj = 0, cex = 1.2, font = 2, col = COL$text)
  mtext(sub,  side = 3, line = 1.2, adj = 0, cex = 0.85, col = COL$text2)
}

# show_ewma = TRUE also draws the EWMA normal VaR and its breaches.
plot_returns_var <- function(fc, file, from = min(fc$date), to = max(fc$date),
                             period_label, date_format = "%Y", by = "year",
                             show_ewma = FALSE) {
  d <- fc[fc$date >= from & fc$date <= to, ]
  r  <- 100 * d$ret_port
  vh <- -100 * d$var_hist
  vn <- -100 * d$var_normal
  ve <- -100 * d$var_ewma_normal
  bh <- d$breach_historical
  bn <- d$breach_normal
  be <- d$breach_ewma_normal

  open_png(file, bottom = if (show_ewma) 7.6 else 6.2)
  plot(d$date, r, type = "n", xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "Daily portfolio return",
       ylim = range(c(r, vh, vn, if (show_ewma) ve)) * 1.04)
  pct_axis()
  ticks <- seq(as.Date(cut(from, by)), to, by = by)
  axis(1, at = ticks, labels = format(ticks, date_format), lwd = 0, lwd.ticks = 1,
       col.ticks = COL$axis)
  abline(h = 0, col = COL$axis, lwd = 1)

  lines(d$date, r, type = "h", col = COL$returns, lwd = 0.7)
  lines(d$date, vn, type = "s", col = COL$normal, lwd = 1.8)
  lines(d$date, vh, type = "s", col = COL$hist, lwd = 1.8)
  if (show_ewma) lines(d$date, ve, type = "s", col = COL$ewma, lwd = 2.2)
  points(d$date[bn], r[bn], pch = PCH$normal, col = COL$normal, cex = 1.5, lwd = 1.8)
  points(d$date[bh], r[bh], pch = PCH$hist, col = COL$hist, cex = 0.95)
  if (show_ewma) {
    # EWMA breach markers sit just left of the return so all three stay visible.
    off <- 0.45 * as.numeric(diff(range(d$date))) / 100
    points(d$date[be] - off, r[be], pch = PCH$ewma, col = COL$ewma, cex = 1.3)
  }

  title_block(
    sprintf("60/40 SPY/AGG portfolio: daily returns vs 1-day 99%% VaR, %s", period_label),
    "VaR is drawn as a negative return threshold. A breach is a day whose loss went below its model's line."
  )
  if (!show_ewma) {
    bottom_legend(c("Daily return",
                    "Historical VaR (250-day window)",
                    sprintf("Historical breach (%d)", sum(bh)),
                    "Normal VaR (250-day window)",
                    sprintf("Normal breach (%d)", sum(bn))),
           cex = 0.85, seg.len = 1.6,
           col = c(COL$returns, COL$hist, COL$hist, COL$normal, COL$normal),
           lty = c(1, 1, NA, 1, NA), lwd = c(2, 2, NA, 2, 1.8),
           pch = c(NA, NA, PCH$hist, NA, PCH$normal))
  } else {
    # Two legend rows: the lines, then the breach markers.
    bottom_legend(c("Daily return", "Historical VaR (250-day window)",
                    "Normal VaR (250-day window)", "EWMA normal VaR (lambda 0.94)"),
                  cex = 0.85, inset = -0.19, seg.len = 1.6,
                  col = c(COL$returns, COL$hist, COL$normal, COL$ewma),
                  lty = 1, lwd = c(2, 2, 2, 2.2))
    bottom_legend(c(sprintf("Historical breach (%d)", sum(bh)),
                    sprintf("Normal breach (%d)", sum(bn)),
                    sprintf("EWMA normal breach (%d)", sum(be))),
                  cex = 0.85, inset = -0.25,
                  col = c(COL$hist, COL$normal, COL$ewma),
                  pch = c(PCH$hist, PCH$normal, PCH$ewma),
                  pt.cex = c(0.95, 1.3, 1.2), pt.lwd = c(1, 1.8, 1), x.intersp = 1)
  }
  dev.off()
}

plot_breaches_by_year <- function(by_year, file = "breaches_per_year.png") {
  models <- names(MODEL_COL)
  m <- t(as.matrix(by_year[, paste0(models, "_breaches")]))
  dimnames(m) <- list(models, by_year$year)
  cols <- unname(MODEL_COL[models])
  space <- c(0.12, 1)

  open_png(file, width = 2800, height = 1300)
  ymax <- max(m, by_year$expected) + 1.5
  bp <- barplot(m, beside = TRUE, col = cols, border = COL$bg,
                space = space, ylim = c(0, ymax), axes = FALSE,
                ylab = "Number of breaches", xlab = "", names.arg = colnames(m),
                col.axis = COL$text2, cex.names = 0.95)
  at <- seq(0, ceiling(ymax), by = 2)
  abline(h = at, col = COL$grid, lwd = 0.8)
  axis(2, at = at, lwd = 0)
  # redraw bars over the grid
  barplot(m, beside = TRUE, col = cols, border = COL$bg,
          space = space, axes = FALSE, names.arg = rep("", ncol(m)), add = TRUE)
  # Expected breaches per year = 1% of that year's trading days (2.50 to 2.53).
  half <- 0.5 * (bp[nrow(bp), 1] - bp[1, 1]) + 0.6
  mid  <- colMeans(bp)
  segments(mid - half, by_year$expected, mid + half, by_year$expected,
           col = COL$text, lwd = 2, lty = 2)

  # Count labels above each bar, on a small background patch so the dashed
  # line never runs through a number.
  ly <- m + 0.45
  w  <- strwidth("00", cex = 0.7) * 0.6
  h  <- strheight("0", cex = 0.7) * 0.85
  rect(bp - w, ly - h, bp + w, ly + h, col = COL$bg, border = NA)
  text(bp, ly, labels = m, cex = 0.7, col = COL$text)

  title_block(
    "VaR breaches per calendar year, 60/40 SPY/AGG portfolio, 2015-2025",
    "Bars: breaches of each 1-day 99% VaR model. Dashed line: expected breaches if the model is correct (1% of trading days, about 2.5 a year)."
  )
  bottom_legend(c(MODEL_LABEL[models], "Expected at 1% of days"), cex = 0.85,
         col = c(cols, COL$text), pch = c(rep(15, length(models)), NA),
         pt.cex = 1.8, lty = c(rep(NA, length(models)), 2),
         lwd = c(rep(NA, length(models)), 2), seg.len = 2)
  dev.off()
}

run_charts <- function(fc, by_year) {
  plot_returns_var(fc, "returns_vs_var_full.png",
                   period_label = "2015-2025")
  plot_returns_var(fc, "returns_vs_var_2020_2021.png",
                   from = as.Date("2020-01-01"), to = as.Date("2021-04-30"),
                   period_label = "January 2020 to April 2021",
                   date_format = "%b %Y", by = "month", show_ewma = TRUE)
  plot_breaches_by_year(by_year)
  message("charts saved to ", FIG_DIR, "/")
}
