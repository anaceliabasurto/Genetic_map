### CONSTRUCCION DEL MAPA GENETICO ###


cross1<-read.cross(format="csv", file = "~/tgbs_data/Genotype_corrector/CR_60/60_2/60_OK/F2_ABHcorrected_final_sin_parentales.csv",
                   estimate.map = FALSE)
summary(cross1)

feno <- read.csv("~/tgbs_data/Genotype_corrector/CR_60/60_2/60_OK/Fenotipos.csv", header = T)
head(feno)




cross1<-convert2bcsft(cross1,BC.gen = 0,F.gen = 0,estimate.map = FALSE)

summary(cross1)


cross1 <- pullCross(cross1, type = "co.located")

supermap<-mstmap(cross1, id = 'id',bychr = TRUE, dist.fun = 'kosambi', p=2)
supermap<-flip.order(supermap,c(1:12))

View(supermap)

write.cross(supermap, "csv", filestem = "file", (1:12))


map<-pull.map(supermap)
plot.map(map)

View(supermap)
summary(supermap)
jittermap(supermap)
plotMissing(supermap)


plotRF(supermap)


supermap <- calc.genoprob(supermap)