setwd("C:/Users/yanni/Desktop/CNB/d1Salmonella")
dfcutoff <- read_excel("cutoffs_def0.001.xlsx")
View(dfcutoff)  

cols <- c("ΔlppA early", "ΔlppA late", "ΔlppB early", "ΔlppB late", "WT early", "WT late")

dfcutoff[cols] <- lapply(dfcutoff[cols], function(x) as.numeric(x))

#ΔlppA early
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `ΔlppA early`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `ΔlppA early`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="ΔlppA early")+
  theme_bw()

#ΔlppA late
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `ΔlppA late`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `ΔlppA late`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="ΔlppA late")+
  theme_bw()

#ΔlppB early
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `ΔlppB early`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `ΔlppB early`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="ΔlppB early")+
  theme_bw()

#ΔlppB late
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `ΔlppB late`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `ΔlppB late`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="ΔlppB late")+
  theme_bw()

#WT early
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `WT early`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `WT early`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="WT early")+
  theme_bw()

#WT late
ggplot(dfcutoff) +
  geom_line(aes(x = Cutoff, y = `WT late`, group = 1)) +
  geom_point(aes(x = Cutoff, y = `WT late`)) +
  labs(y="Valor", 
       x="Cutoffs",
       title="WT late")+
  theme_bw()
