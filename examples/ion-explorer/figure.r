library(ggplot2)

for (file in c("start", "fragment_b", "fragment_y", "fragment_v_0", "fragment_w_0", "fragment_w_1", "fragment_w_2", "fragment_w_3", "fragment_w_4", "fragment_w_5", "fragment_w_6")) {
    data = read.csv(paste(file, ".csv", sep=""), header=TRUE);
    # data = data[data$count > 1000,];

    plot = ggplot(data, aes(x=mz, y=count, size=avg_intensity, colour=avg_intensity)) +
        geom_point() +
        scale_size_continuous(range = c(0.25, 3)) +
        theme_bw() 

    if (file == "fragment_y") {
        plot = plot + geom_vline(xintercept=0, linetype="dashed") + geom_vline(xintercept=-18.011+0.984, linetype="dashed") + geom_vline(xintercept=+25.979-18.011, linetype="dashed")
    }

    ggsave(paste(file, ".png", sep=""), plot)
}

for (file in c("comparison_b_N_N[U:Deamidated]", "comparison_y_N_N[U:Deamidated]")) {
    data = read.csv(paste(file, ".csv", sep=""), header=TRUE);
    data = data[data$count_a > 20 | data$count_b > 20,];
    max_a = max(data$avg_intensity_a);
    max_b = max(data$avg_intensity_b);

    ggplot(data, aes(x=mz, y=log2(count_a/count_b), colour=rgb(avg_intensity_a / max_a, 0, avg_intensity_b / max_b), size=avg_intensity_a)) +
        geom_point() +
        scale_size_continuous(range = c(0.25, 3)) +
        theme_bw() +
        guides(colour="none")

    ggsave(paste(file, ".png", sep=""))
}
