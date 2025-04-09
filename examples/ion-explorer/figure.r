library(ggplot2)
library(dplyr)

# for (file in c("start", "fragment_b", "fragment_y", "fragment_v_0", "fragment_w_0", "fragment_w_1", "fragment_w_2", "fragment_w_3", "fragment_w_4", "fragment_w_5", "fragment_w_6", "fragment_precursor")) {
#     data = read.csv(paste("examples/ion-explorer/data/", file, ".csv", sep=""), header=TRUE);
#     # data = data[data$count > 1000,];

#     plot = ggplot(data, aes(x=mass, y=count, size=avg_intensity, colour=avg_intensity)) +
#         geom_point() +
#         scale_size_continuous(range = c(0.25, 3)) +
#         xlab("Difference to theoretical ion (Da)") +
#         theme_bw() 

#     if (file == "fragment_y") {
#         plot = plot + geom_vline(xintercept=0, linetype="dashed") + geom_vline(xintercept=+18.011-0.984, linetype="dashed") + geom_vline(xintercept=-25.979+18.011, linetype="dashed")
#     }

#     ggsave(paste("examples/ion-explorer/figures/", file, ".png", sep=""), plot)
# }

merge_stack = function(a, b) {
    # print(head(a));
    # print(head(b));

    return(rbind(a, b) %>%
        group_by(mass) %>%
        reframe(
            count = sum(count, na.rm = TRUE),
            total_intensity = sum(total_intensity, na.rm = TRUE),
        ));
}

get_stack = function(fragment, element=NA, mode=NA) {
    pattern = paste("examples/ion-explorer/data/fragment_", (if (is.na(fragment)) "*" else fragment), "_", (if (is.na(element)) "*" else element), "_", (if (is.na(mode)) "*" else mode), ".csv", sep="");
    files = Sys.glob(pattern);
    # print(files);
    stack = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("mass", "count", "total_intensity"));
    for (file in files) {
        data = read.csv(file, header=TRUE);
        stack = merge_stack(stack, data);
    }
    return(stack);
}

make_plot = function(fragment, element=NA, mode=NA) {
    data = get_stack(fragment, element, mode);
    plot = ggplot(data, aes(x=mass, y=count, size=total_intensity / count, colour=total_intensity / count)) +
        geom_point() +
        scale_size_continuous(range = c(0.25, 3)) +
        xlab("Difference to theoretical ion (Da)") +
        theme_bw() 

    if (fragment == "y") {
        plot = plot + geom_vline(xintercept=0, linetype="dashed", alpha=0.25) + geom_vline(xintercept=0.984-18.011, linetype="dashed", alpha=0.25) + geom_vline(xintercept=25.979-18.011, linetype="dashed", alpha=0.25) + xlim(-50, 50)
    }
    if (fragment == "b") {
        plot = plot + geom_vline(xintercept=0, linetype="dashed", alpha=0.25) + geom_vline(xintercept=-27.995, linetype="dashed", alpha=0.25) + geom_vline(xintercept=+17.027, linetype="dashed", alpha=0.25) + xlim(-50, 50)
    }

    ggsave(paste("examples/ion-explorer/figures/fragment_", (if (is.na(fragment)) "-" else fragment), "_", (if (is.na(element)) "-" else element), "_", (if (is.na(mode)) "-" else mode), ".png", sep=""), plot)
}

compare_plot = function(fragment1, fragment2, element1, element2, mode1, mode2) {
    data1 = get_stack(fragment1, element1, mode1);
    data2 = get_stack(fragment2, element2, mode2);
    max1 = max(data1$count);
    max2 = max(data2$count);
    # merged = merge(by = "mass", all = TRUE, x = data1, y = data2, no.dups = TRUE);
    # print(head(data1));
    # print(head(data2));

    plot = ggplot(data1) +
        geom_point(aes(x=mass, y=count / max1, size=total_intensity / count), colour="blue") +
        geom_point(data=data2, aes(x=mass, y=-count / max2, size=total_intensity / count), colour="red") +
        scale_size_continuous(range = c(0.125, 3)) +
        xlab("Difference to theoretical ion (Da)") +
        theme_bw() 

    if (fragment1 == "y" && fragment2 == "y") {
        plot = plot + 
            geom_vline(xintercept=0, linetype="dashed", alpha=0.25) + 
            geom_vline(xintercept=0.984-18.011, linetype="dashed", alpha=0.25) + 
            geom_vline(xintercept=0.984-18.011-15.023, linetype="dotted", alpha=0.25) + 
            geom_vline(xintercept=0.984-18.011-29.039, linetype="dotted", alpha=0.25) + 
            geom_vline(xintercept=0.984-18.011-43.088, linetype="dotted", alpha=0.25) + 
            geom_vline(xintercept=25.979-18.011, linetype="dashed", alpha=0.25) + 
            xlim(-75, 50)
    }
    if (fragment1 == "b" && fragment2 == "b") {
        plot = plot + geom_vline(xintercept=0, linetype="dashed", alpha=0.25) + geom_vline(xintercept=-27.995, linetype="dashed", alpha=0.25) + geom_vline(xintercept=+17.027, linetype="dashed", alpha=0.25) + xlim(-75, 50)
    }

    ggsave(paste("examples/ion-explorer/figures/compare_", (if (is.na(fragment1)) "-" else fragment1), "_", (if (is.na(element1)) "-" else element1), "_", (if (is.na(mode1)) "-" else mode1), "_to_", (if (is.na(fragment2)) "-" else fragment2), "_", (if (is.na(element2)) "-" else element2), "_", (if (is.na(mode2)) "-" else mode2), ".png", sep=""), plot, unit="cm", width= 32, height=18)
}

# for (fragment in c("y", "b", "precursor")) {
#     for (mode in c("EAciD", "EAD", "CID")) {
#         make_plot(fragment, NA, mode);
#     }
# }

compare_plot("y", "y", "D", "N(Deamidated)", NA, NA);
compare_plot("y", "y", "I", "L", "EAciD", "EAciD");
compare_plot("b", "b", "D", "N(Deamidated)", NA, NA);
compare_plot("b", "b", "I", "L", NA, NA);