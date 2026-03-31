# ==============================================================================
# Gráfico — ARI, AMI, NMI por dataset (un panel por dataset)
# ==============================================================================

resumen_dataset <- datos |>
  group_by(algoritmo, name) |>
  summarise(
    ARI = mean(ARI, na.rm = TRUE),
    AMI = mean(AMI, na.rm = TRUE),
    NMI = mean(NMI, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(ARI, AMI, NMI),
               names_to  = "metrica",
               values_to = "valor")

ggplot(resumen_dataset, aes(x = reorder(algoritmo, valor), 
                            y = valor, 
                            fill = metrica)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = round(valor, 2)),
    position = position_dodge(width = 0.8),
    hjust    = -0.15,
    size     = 3.5        # ← más grande
  ) +
  coord_flip() +
  facet_wrap(~ name, ncol = 1) +     # ← un panel por dataset
  scale_fill_manual(values = c(ARI = "#4C72B0", AMI = "#DD8452", NMI = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Métricas externas por algoritmo y dataset",
    x        = "Algoritmo",
    y        = "Valor",
    fill     = "Métrica"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "top",
    strip.text       = element_text(face = "bold", size = 11),  # título de cada panel
    strip.background = element_rect(fill = "gray92", color = NA)
  )

ggsave("grafico_metricas_por_dataset.png", width = 12, height = 7 * n_distinct(datos$name), 
       dpi = 300, limitsize = FALSE)

# ==============================================================================
# Gráfico — Silhouette por dataset (un panel por dataset)
# ==============================================================================

resumen_sil_dataset <- datos |>
  group_by(algoritmo, name) |>
  summarise(
    Silhouette = mean(Silhouette_mean, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(resumen_sil_dataset, aes(x = reorder(algoritmo, Silhouette),
                                y = Silhouette,
                                fill = Silhouette > 0)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray40", linewidth = 0.5) +
  geom_text(
    aes(label = round(Silhouette, 3),
        hjust = ifelse(Silhouette >= 0, -0.15, 1.15)),
    size = 3
  ) +
  coord_flip() +
  facet_wrap(~ name, ncol = 1) +     # ← un panel por dataset
  scale_fill_manual(values = c("TRUE" = "#55A868", "FALSE" = "#C44E52"),
                    labels = c("TRUE" = "Positivo", "FALSE" = "Negativo")) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  labs(
    title    = "Coeficiente de Silhouette por algoritmo y dataset",
    subtitle = "Verde = clusters bien formados  |  Rojo = clusters mal formados",
    x        = "Algoritmo",
    y        = "Silhouette promedio",
    fill     = "Valor"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "top",
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray92", color = NA)
  )

ggsave("grafico_silhouette_por_dataset.png", width = 12, 
       height = 5 * n_distinct(datos$name),
       dpi = 300, limitsize = FALSE)