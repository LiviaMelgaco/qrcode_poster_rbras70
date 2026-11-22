##### GRAFICOS POSTER ######

library(ggplot2)
library(scales)

dados.plot <- data.frame(
  phi    = rep(c("φ = 0,3", "φ = 0,8"), times = 2),
  metodo = rep(c("GLARMA clássico", "GLARMA bayesiano"), each = 2),
  eqm    = c(0.0069, 0.1724,   # clássico: phi = 0,3 e 0,8
             0.0071, 0.0667)   # bayesiano: phi = 0,3 e 0,8
)
dados.plot$metodo <- factor(dados.plot$metodo,
                       levels = c("GLARMA clássico", "GLARMA bayesiano"))

g <- ggplot(dados.plot, aes(x = phi, y = eqm, fill = metodo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = format(eqm, decimal.mark = ",")),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("GLARMA clássico"  = "#9AA5B1",
                               "GLARMA bayesiano" = "#1F4E8C")) +
  scale_y_continuous(limits = c(0, 0.2),
                     labels = label_number(decimal.mark = ",")) +
  labs(x = NULL, y = "EQM de beta 0", fill = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank())

g
ggsave("fig1_eqm_beta0.png", g, width = 4, height = 3, dpi = 300)
