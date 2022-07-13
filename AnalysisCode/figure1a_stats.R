print(Sys.time())

setwd("/mBE/")

marks <- c("H2AZ", "H2BK120ac", "H2BK12ac", "H2BK5ac", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3", "H4K20me1", "H4K8ac", "mH2A1", "mH2A2")
HMEC_bed <- read.table("/mBE/E119_25_imputed12marks.bed", sep="\t")
HMEC1 <- read.table(file = "/mBE/HMEC_25st_alt.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
HMEC1 <- HMEC1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
NHM_bed <- read.table("/mBE/E059_25_imputed12marks.bed", sep="\t")
NHM1 <- read.table(file = "/mBE/NHM_25st.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
NHM1 <- NHM1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HepG2_bed <- read.table("/mBE/E118_25_imputed12marks.bed", sep="\t")
HepG2_1 <- read.table(file = "/mBE/HepG2_25st.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
HepG2_1 <- HepG2_1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HMEC_bed <- HMEC_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HMEC_states <- HMEC_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)
NHM_bed <- NHM_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
NHM_states <- NHM_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)
HepG2_bed <- HepG2_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HepG2_states <- HepG2_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)

combos <- expand.grid(as.character(1:25), marks)
combos_split <- split(combos, seq(nrow(combos)))
get_tests <- function(combo){
    x <- curr %>% filter(state == as.integer(combo$Var1)) %>% select(combo$Var2) %>% unlist() %>% unname()
    y <- curr %>% filter(state != as.integer(combo$Var1)) %>% select(combo$Var2) %>% unlist() %>% unname()
    p <- wilcox.test(x = x, y = y)$p.value
    padj <- p * 400 # Bonferroni correction
    return(padj >= 0.05)
}
curr <- HMEC1
HMEC_insignif <- unname(unlist(lapply(combos_split, get_tests)))
curr <- NHM1
NHM_insignif <- unname(unlist(lapply(combos_split, get_tests)))
curr <- HepG2_1
HepG2_insignif <- unname(unlist(lapply(combos_split, get_tests)))

insignif <- bind_rows(list(HMEC = combos[HMEC_insignif,], NHM = combos[NHM_insignif,], HepG2 = combos[HepG2_insignif,]), .id = "cell")
write.table(insignif, file = "fig1a_stats.tab", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

print(Sys.time())
