require(ggplot2)
suppressMessages(require(magrittr))
require(cowplot)
require(MASS)

set.seed(1212111)
sigma <- rbind(c(1, 0.7), c(0.7, 1))
mu <- c(10, 5)


x <- as.data.frame(mvrnorm(n = 200, mu = mu, Sigma = sigma))

distx <- x[,1] + rnorm(200)
## hist(x[,1])
## hist(distx)

alldata <- rbind(x, data.frame(V1 = distx, V2 = x[,2]))
alldata$type <- as.factor(c(rep("No error", 200), rep("Error", 200)))

alldata %>% ggplot(aes(x = V1, fill = type, group = type)) + geom_density(adjust = 1.5, alpha=.4) + theme_minimal() + xlab("Y") + scale_fill_brewer(palette="Dark2") -> spread

alldata %>% ggplot(aes(x = V2, y = V1, col = type, group = type)) + geom_point(alpha = .4) + geom_smooth(aes(group=type), method = "lm", formula = 'y ~ x', alpha = .4, se = FALSE) + theme_minimal() + xlab("X") + ylab("Y") + scale_color_brewer(palette="Dark2") -> scatter

plot_grid(spread, scatter)
