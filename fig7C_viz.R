library(ggplot2)
library(ggpubr)

#Experimental replicate data: Emax from Fig7B, p1111 from Fig7A
corr_data <- data.frame(
  UV_min = factor(rep(c(0,15,30,60), times=3)),
  Emax = c(0.005566,-0.3221,-0.4886,-0.7124,
           0.06547,-0.2288,-0.2991,-1.178,
           0.02376,-0.1757,-0.3198,-0.9587),
  p1111 = c(0.4725878,0.3899299,0.3709433,0.3239317,
            0.4681347,0.4119134,0.4037808,0.289764,
            0.407683,0.3786112,0.3649559,0.240711)
)

#Pearson correlation
cor_test <- cor.test(corr_data$p1111, corr_data$Emax)
cor_test

Fig7C <- ggplot(corr_data, aes(x=p1111, y=Emax)) +
  geom_point(aes(color=UV_min), size=4) +
  scale_color_manual(values=c("0"="#A0A0A4", "15"="#1CC5FE", "30"="#C000C0", "60"="#FF8000"))+
  geom_smooth(method="lm", se=TRUE, color="black", fill="lightgrey", alpha=0.3)+
  scale_x_continuous(limits = c(0.2, 0.5))+
  scale_y_continuous(limits = c(-1.5, NA))+
  labs(x="Probability Intact",
       y="Emax",
       color="UV Time (min)") +
  theme_classic2(base_size = 12)

lm_model <- lm(Emax ~ p1111, data = corr_data)
summary(lm_model)$adj.r.squared


ggsave(plot = Fig7C, filename="fig7_corr_panel.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4)

ggsave(plot = Fig7C, filename="fig7_corr_panel.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4, 
       dpi=600,
       bg="white")




#LMM replicate effect comparison
library(lme4)      
library(lmerTest)  

lmm_data <- data.frame(
  replicate = factor(rep(1:3, each=4)),
  UV_min = factor(rep(c(0,15,30,60), times=3)),
  Emax = c(0.005566,-0.3221,-0.4886,-0.7124,
           0.06547,-0.2288,-0.2991,-1.178,
           0.02376,-0.1757,-0.3198,-0.9587),
  p1111 = c(0.4725878,0.3899299,0.3709433,0.3239317,
            0.4681347,0.4119134,0.4037808,0.289764,
            0.407683,0.3786112,0.3649559,0.240711)
)

lmm <- lmer(Emax ~ p1111 + (1|replicate), data=lmm_data)
summary(lmm)  

#compare both

# Predicted values for LMM fixed effect
pred_data <- data.frame(p1111 = seq(min(lmm_data$p1111), max(lmm_data$p1111), length.out=100))
pred <- predict(lmm, newdata=pred_data, re.form=NA, se.fit=TRUE)
pred_data$fit <- pred$fit
pred_data$lower <- pred$fit - 1.96*pred$se.fit
pred_data$upper <- pred$fit + 1.96*pred$se.fit

# LM predictions
lm_model <- lm(Emax ~ p1111, data=lmm_data)
lm_pred <- data.frame(p1111 = pred_data$p1111,
                      fit = predict(lm_model, newdata=pred_data))

ggplot(lmm_data, aes(x=p1111, y=Emax, color=UV_min)) +
  geom_point(size=4) +
  scale_color_manual(values=c("0"="#A0A0A4", "15"="#1CC5FE",
                              "30"="#C000C0", "60"="#FF8000")) +
  geom_line(data=lm_pred, aes(x=p1111, y=fit), color="blue", size=1, linetype="dashed", inherit.aes = FALSE) +
  geom_line(data=pred_data, aes(x=p1111, y=fit), color="black", size=1, inherit.aes = FALSE) +
  theme_classic2(base_size=12)

