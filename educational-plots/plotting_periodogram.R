library(ggplot2)
library(reshape2)  # To transform the matrix for plotting
library(latex2exp) # For LaTeX math rendering in titles

# --- 1. Define the Plotting Function ---
plot_academic_matrices <- function(mat_original, mat_shuffled, title1, title2) {
  # --- 1. Coordinate Calculation (Fourier Frequencies) ---
  # Calculate N for rows and cols (assuming square, but this works for rect too)
  N_row <- nrow(mat_original)
  N_col <- ncol(mat_original)

  # Create the sequence from -floor((N-1)/2) to floor(N/2)
  freqs_row <- seq(-floor((N_row - 1)/2), floor(N_row/2))
  freqs_col <- seq(-floor((N_col - 1)/2), floor(N_col/2))

  df_orig <- melt(as.matrix(mat_original))

  df_orig$Label <- TeX(title1)

  df_shuf <- melt(as.matrix(mat_shuffled))
  df_shuf$Label <- TeX(title2)

  df_combined <- rbind(df_orig, df_shuf)

  # Factor to ensure correct ordering
  df_combined$Label <- factor(df_combined$Label,
                              levels = c(TeX(title1),
                                         TeX(title2)))
  # --- 3. Map Indices to Frequencies ---
  df_combined$Freq_u <- freqs_row[df_combined$Var1]
  df_combined$Freq_v <- freqs_col[df_combined$Var2]

  p <- ggplot(df_combined, aes(x = Freq_u, y = Freq_v, fill = value)) +
    geom_raster() +
    coord_fixed() +

    # label_parsed will now understand the strings because we added quotes and ~
    facet_wrap(~ Label, ncol = 2, labeller = label_parsed) +

    scale_fill_distiller(direction = 1, name = "Amplitude") +

    theme_minimal(base_size = 16, base_family = "Times") +

    labs(
      x = TeX("$\\omega_l$"),
      y = TeX("$\\omega_k$")
    ) +

    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 25, hjust = 0.5, margin = margin(b=15)),
      axis.title = element_text(size = 20),
      panel.grid = element_blank(),
      legend.position = "right",
      strip.text = element_text(size = 25)
    )

  print(p)
}

x <- gridMA(100, 100, MA_coef_row(.7))
y <- gridMA(100, 100, MA_coef_col(.7))

I_x <- I(x)
I_y <- I(y)

I_tilde <- .5 * (I_x + I_y)
I_x_diff <- I_x - I_tilde
I_y_diff <- I_y - I_tilde

perm <- generate_random_mask(100, 100)
perm_n <- matrix(as.numeric(perm != 1), 100, 100)

I_x_rand <- I_x_diff * perm + I_y_diff * perm_n
I_y_rand <- I_y_diff * perm + I_x_diff * perm_n
h <- (100^2)^(-.2)
Kh <- k_2d_bp(100, 100, h, h)

plot_academic_matrices(I_x_diff,
                       I_y_diff,
                       "$I_x -\\tilde{I}$",
                       "$I_y - \\tilde{I}$")

plot_academic_matrices(I_x_rand,
                       I_y_rand,
                       "$I_x^\\pi -\\tilde{I}$",
                       "$I_y^\\pi - \\tilde{I}$")



plot_academic_matrices(
  EBImage::filter2(
    I_x_diff, Kh, "circular"
  )^2,
  EBImage::filter2(
    I_y_diff, Kh, "circular"
  )^2,
  "$(\\hat{f}_x -\\tilde{f})^2$",
  "$(\\hat{f}_y - \\tilde{f})^2$")


plot_academic_matrices(
  EBImage::filter2(
    I_x_rand, Kh, "circular"
  )^2,
  EBImage::filter2(
    I_y_rand, Kh, "circular"
  )^2,
  "$(\\hat{f}_x^\\pi -\\tilde{f})^2$",
  "$(\\hat{f}_y^\\pi - \\tilde{f})^2$")
