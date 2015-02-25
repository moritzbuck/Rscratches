
pdf("bin_125.pdf", width=8, height=8)
plot(plot.cov.pcas(fact=1, bins=c("125"), c.list=data$clusts %in% bin.stats$clusts[bin.stats$ass_size > 0.5] )+scale_y_log10())
dev.off()

pdf("bin_125_time_series.pdf", width=10, height=4)
plot(plot.tseries(bins=c("125")))
dev.off()
