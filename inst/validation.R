library(lidRalignment)
library(dplyr)
library(tidyr)

files = list.files("/home/jr/Documents/Usherbrooke/data/ForestReg", recursive = TRUE, pattern = ".las", full.names = TRUE)
dirs <- dirname(files)
plot <- basename(dirs)
difficulty <- basename(dirname(dirs))
files <- data.frame(
  difficulty = difficulty,
  plot = plot,
  filename = files,
  stringsAsFactors = FALSE
)
files$type <- ifelse(startsWith(basename(files$filename), "TLS"), "mov", ifelse(startsWith(basename(files$filename), "ULS"), "fix", NA))
files <- pivot_wider(files, id_cols = c("difficulty", "plot"), names_from = "type", values_from = "filename")

matrices = list.files("/home/jr/Documents/Usherbrooke/data/ForestReg", recursive = TRUE, pattern = ".txt", full.names = TRUE)
dirs <- dirname(matrices)
plot <- basename(dirs)
difficulty <- basename(dirname(dirs))
matrices <- data.frame(
  difficulty = difficulty,
  plot = plot,
  mat = matrices,
  stringsAsFactors = FALSE
)
matrices$type <- ifelse(startsWith(basename(matrices$mat), "ground"), "Mt", ifelse(startsWith(basename(matrices$mat), "matrice"), "M", NA))
matrices <- pivot_wider(matrices, id_cols = c("difficulty", "plot"), names_from = "type", values_from = "mat")

data = left_join(files, matrices, by = c("difficulty", "plot"))

# i = 2 FAILURE
# S01 i = 1
# S02 i = 2
# S04 i = 3
# S05 i = 12  si R trop grand
# S06 i = 4 OK après ICP 3 axes
# S07 i = 5
#i = 14 bug
# S08 i = 13  21 m d'erreur BUG
# S10 i = 6
# S14 i = 9
# S15 i = 17
i = 1

for (i in 1:17)
{
  fix = data$fix[i]
  mov =  data$mov[i]
  mat =  data$Mt[i]
  Mt = as.matrix(data.table::fread(mat))
  Mt

  err = max(abs(Mt[1:2,4]))
  max_offset = 10
  if (err > 10) max_offset = err + 2


  alignment = AlignmentScene$new(fix, mov)
  alignment$set_ref_is_ground_based(FALSE)
  alignment$set_mov_is_ground_based(TRUE)
  alignment$set_radius(30)

  alignment$prepare()
  #alignment$plot("raw")

  alignment$coarse_align(res = 2, max_offset = max_offset, debug = F)
  #alignment$plot("coarse")

  alignment$fine_align(use_cc = F)
  #alignment$plot("fine")

  M = alignment$get_registration_matrix()
  M

  omat = paste0(dirname(mat), "/matrice.txt")
  data.table::fwrite(as.data.frame(M), omat)
}

compare_matrices = function(M, Mt)
{
  dx = M[1,4] - Mt[1,4]
  dy = M[2,4] - Mt[2,4]
  dz = M[3,4] - Mt[3,4]
  rzmrad  = (atan2(M[2,1], M[1,1]) - atan2(Mt[2,1], Mt[1,1])) * 1000
  rzdeg = rzmrad/1000*180/pi
  dxy = sqrt(dx^2+dy^2)*100
  dz = sqrt(dz^2)*100
  cat("dxy =", round(dxy,1) , "cm\ndz =", round(dz,1), "cm\nrz =", round(rzmrad, 1), "mrad\nrz =", round(rzdeg, 1), "deg\n")
  return(list(dxy = dxy, dz = dz, rzmrad = rzmrad, rzdeg = rzdeg))
}

u = compare_matrices(M, Mt)

if (FALSE)
{
  ref = lidR::readLAS(fix, select = "xyz", filter = "-keep_random_fraction 0.2")
  mov = lidR::readLAS(mov, select = "xyz", filter = "-keep_random_fraction 0.1")
  ref@header$`X scale factor` = 0.001
  ref@header$`Y scale factor` = 0.001
  ref@header$`Z scale factor` = 0.001
  lidR::las_quantize(ref)
  mov@header$`X scale factor` = 0.001
  mov@header$`Y scale factor` = 0.001
  mov@header$`Z scale factor` = 0.001
  lidR::las_quantize(mov)
  lidRalignment:::show_alignment(ref, mov, Mt)
}


ans = apply(data, 1, function(x)
{
  Mt = as.matrix(data.table::fread(x[5]))
  M = as.matrix(data.table::fread(x[6]))
  compare_matrices(M, Mt)
})

ans = data.table::rbindlist(ans)
res = cbind(data, ans)

library(ggplot2)
library(gridExtra)
library(plotly)

fail = filter(res, abs(rzdeg) > 10)
res = filter(res, abs(rzdeg) < 10)

nolegend = theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 <- ggplot(res) +
  aes(x = plot, y = dxy, col = difficulty) +
  geom_point() +
  xlab("Plot") + ylab("Avegrage translation residual (cm)")+
  nolegend
p2 <- ggplot(res) +
  aes(x = plot, y = dz, col = difficulty) +
  geom_point() +
  nolegend
p3 <- ggplot(res) +
  aes(x = plot, y = abs(rzdeg), col = difficulty) +
  geom_point() +
  xlab("Plot") + ylab("Avegrage rotation residual (deg)") +
  nolegend

grid.arrange(p1, p3, ncol = 2)  # side by side
#grid.arrange(p1, p2, p3, ncol = 2)  # side by side

plotly::ggplotly(p1)

p1 <- ggplot(res) +
  aes(x = plot, y = dxy, group = difficulty, fill = difficulty) +
  geom_boxplot() +
  geom_hline(yintercept = c(50, 100), lty = 3) +
  xlab("Plot") + ylab("Avegrage translation residual (cm)")
p1
p2 <- ggplot(res) +
  aes(x = plot, y = abs(rzdeg), group = difficulty, fill = difficulty) +
  geom_boxplot() +
  geom_hline(yintercept = c(1, 2), lty = 3) +
  xlab("Plot") + ylab("Avegrage rotation residual (deg)")
p2

group_by(res, difficulty) %>%  summarise(trans = mean(dxy), rot = mean(abs(rzdeg)))

group_by(res, difficulty) %>%  summarise(trans = median(dxy), rot = median(abs(rzdeg)))
