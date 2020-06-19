#This code is for running for reading in data, cleaning data, setting up data frames 
#for analysis, assessing correlations, running the GLMs and t-tests, doing PERMANOVAs and indicator analysis,
#and doing the reserve simulation 
#from Osgood et al. 2020 AfrMarSci

#see Umbrella.cooccur.R for co-occurrence analysis

library(vegan) # for shannon diversity function
library(labdsv) # for DLI values 

#Read in data
BRUV.MaxN<-read.csv("data/Osgood_etal_2020_MaxN_Data.csv", header=TRUE, stringsAsFactors=FALSE)
env<-read.csv("data/Osgood_etal_2020_BRUV_sampling_metadata_environmentals.csv", header=TRUE, stringsAsFactors=FALSE)

########################################################################
##Remove species (ie invertebrates) that are not relevant/ are rare#####
########################################################################
#Remove krill from analysis
BRUV.MaxN<-BRUV.MaxN[, -which(names(BRUV.MaxN)=="Euphausiids")]

#Remove pelagic larva/bait fish rarely seen, sea pens, and worms
BRUV.MaxN<-BRUV.MaxN[, -which(names(BRUV.MaxN)=="Unknown.fish")]
BRUV.MaxN<-BRUV.MaxN[, -which(names(BRUV.MaxN)=="Larval.fish")]
BRUV.MaxN<-BRUV.MaxN[, -which(names(BRUV.MaxN)=="Bait.fish...sardine")]

#Remove jellyfish, sea urchins, sea anemones, sea cucumbers, sea pens, and worms
BRUV.MaxN<-BRUV.MaxN[, -which(names(BRUV.MaxN)%in%c("Snail", 
  "Feathery.sea.pen", "Fanworm", 
    "Eunicid.worm", "Sea.anemone", "Starfish", "Seastar", 
    "Spiny.starfish", "Sea.cucumber", 
    "Cape.urchin", "Helmet.shell", "Sea.star", "Papery.burnupena", 
    "Sea.slug", "Tick.shell", "Red.sea.star","False.plume.anemone",
    "Topshell.snail", "Topshell", "Turban.shell", "Fat.plough.shell", 
    "Root.mouthed.jellyfish", "Sea.urchin", "Jellyfish", "Purple.lipped.dog.whelk" ))]

############################################################################################
#####Setting data frames of species richness, pres/ab, and mean abundance by site###########
############################################################################################

MaxN<-BRUV.MaxN[,-c(1:5)] # Creates abundance only data frame (remove "meta-data" columns)
#head(MaxN) #Check that I was successful by looking at first few rows

BRUV.Richness<-apply(MaxN, 2, function(x) {as.numeric(x>0)}) # Create data frame that is presence/absence (0 and 1) only

#make dataset of MaxN based on site rather than BRUV drop (with MaxN averaged per site)

MaxN.site<-aggregate(MaxN, by=list(env$Site, env$Protection), mean) #mean of MaxN by site per species
MaxN.site<-MaxN.site[,-c(1:2)] # Remove site and protection columns so to make pure site by species matrix

#make dataset of species richness based on site rather than BRUV drop (with MaxN averaged per site)

Site.Richness1<-aggregate(BRUV.Richness, by=list(env$Site, env$Protection), sum) # total number of times each species found per site

Site<-Site.Richness1$Group.1 #saves index of which site applies to each row
Protection<-Site.Richness1$Group.2 #saves index of which protectino level applies to each row

Site.Richness<-apply(Site.Richness1[,-c(1:2)], 2, function(x) {as.numeric(x>0)}) # Returns data frame to 0 and 1 only - ie yes if species EVER found there, no if NEVER found there

#get shark abundance on average by site
shark.abun<-rowSums(MaxN[,which(colnames(MaxN)%in%all.chon)])

Site.SharkAbun<-aggregate(shark.abun, by=list(env$Site, env$Protection), mean)$x # mean total number of chondrichthyans found per site

############################################################################################
#####Setting the names of various ways of breaking up chondrichthyan class #################
############################################################################################
#so can easily subset based on these divisions later - need periods as R deletes the spaces
all.chon<-c("Dark.shyshark", "Puffadder.shyshark", "Broadnose.sevengill.shark", 
  "Leopard.catshark","Pyjama.shark", "Spearnose.skate", "Biscuit.skate", 
                  "Eagle.ray", "Spotted.gully.shark",
                  "Bronze.whaler", "Smooth.hammerhead.shark", 
                  "Smooth.hound.shark",  "St..Joseph.shark", "Tiger.catshark", 
                    "Lesser.guitarfish", "Soupfin.shark", 
                    "Shortnose.spurdog","Short.tail.stingray")

all.shyshark<-c("Dark.shyshark", "Puffadder.shyshark", "Leopard.catshark", 
  "Pyjama.shark", "Tiger.catshark")

all.large.shark<-c("Broadnose.sevengill.shark", 
          "Spotted-gully.shark",
               "Bronze.whaler", "Smooth.hammerhead.shark", 
               "Smooth-hound.shark",  
                "Soupfin.shark")

all.batoid<-c("Lesser.guitarfish", "Eagle.ray", "Biscuit.skate", 
  "Spearnose.skate", "Short.tail.stingray")

sp.of.interest<-c("Red.steenbras", "Galjoen", "Geelbek", "Red.roman", 
  "Red.stumpnose", "Janbruin", "Blacktail", "Shad", "Silver.kob", "Slender.baardman",
  "West.coast.rock.lobster", "White.stumpnose", "Yellowbelly.rock.cod", "Zebra")


#which sites have different divisions of chondrichthyans
#vector of 0 if no large sharks found there ever, 1 if a large shark has been found at that site
large.presab<-as.numeric(rowSums(Site.Richness[,which(colnames(Site.Richness)%in%all.large.shark)])>0)
#vector of 0 if no catsharks found there ever, 1 if a catshark has been found at that site
shyshark.presab<-as.numeric(rowSums(Site.Richness[,which(colnames(Site.Richness)%in%all.shyshark)])>0)
#vector of 0 if no batoids found there ever, 1 if a batoid has been found at that site
batoid.presab<-as.numeric(rowSums(Site.Richness[,which(colnames(Site.Richness)%in%all.batoid)])>0)
#vector of 0 if less than three individual chondrichthyans found there on average, 
#1 if more than three individual chondrichthyans seen at that site on average
manychon.presab<-as.numeric(Site.SharkAbun>3)


#Data frames for assessing the effects of the presence of large sharks, catsharks (shysharks), batois,
#or any chondrichthyan on the abundance and diversity of other taxa
#Therefore these data frames include ALL species BUT the species group whose effect is being analyzed
#ie Max.all.large includes all species BUT large sharks

#data frame averaged by site of MaxN of all species but large sharks:
MaxN.all.large1<-MaxN[,which(!colnames(MaxN)%in%all.large.shark)]
MaxN.all.large<-aggregate(MaxN.all.large1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame average by site of shannon diversity of all species but large sharks:
MaxN.all.large.diver1<-vegan::diversity(MaxN.all.large1, index="shannon")
MaxN.all.large.diver<-aggregate(MaxN.all.large.diver1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame averaged by site of MaxN of all species but catsharks:
MaxN.all.shyshark1<-MaxN[,which(!colnames(MaxN)%in%all.shyshark)]
MaxN.all.shyshark<-aggregate(MaxN.all.shyshark1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame average by site of shannon diversity of all species but catsharks:
MaxN.all.shyshark.diver1<-vegan::diversity(MaxN.all.shyshark1, index="shannon")
MaxN.all.shyshark.diver<-aggregate(MaxN.all.shyshark.diver1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame averaged by site of MaxN of all species but batoids:
MaxN.all.batoid1<-MaxN[,which(!colnames(MaxN)%in%all.batoid)]
MaxN.all.batoid<-aggregate(MaxN.all.batoid1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame average by site of shannon diversity of all species but batoids:
MaxN.all.batoid.diver1<-vegan::diversity(MaxN.all.batoid1, index="shannon")
MaxN.all.batoid.diver<-aggregate(MaxN.all.batoid.diver1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#data frame averaged by site of MaxN of all species but any chondrichthyan:
MaxN.all.chon1<-MaxN[,which(!colnames(MaxN)%in%all.chon)]
MaxN.all.chon<-aggregate(MaxN.all.chon1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

#breaking up MaxN.all.chon by protection level:
MaxN.all.protected<-subset(MaxN.all.chon, Protection=="Protected")
MaxN.all.unprotected<-subset(MaxN.all.chon, Protection=="Unprotected")

#shannon diversity of everything but chondrichthyans by site
diver.all.protected<-vegan::diversity(MaxN.all.protected, index = "shannon") 
diver.all.unprotected<-vegan::diversity(MaxN.all.unprotected, index = "shannon") 
diver.all<-vegan::diversity(MaxN.all.chon, index = "shannon") 

#breaking up chondrichthyan abundance vector by protection level:
shark.abun.protected<-subset(Site.SharkAbun, Protection=="Protected")
shark.abun.unprotected<-subset(Site.SharkAbun, Protection=="Unprotected")

#get total chondrichthyan species richness by site
shark.rich<-rowSums(Site.Richness[,which(colnames(Site.Richness)%in%all.chon)]) #total number of chondrichthyan species at a site

#get total chondrichthyan species richness by site and by protection:
shark.rich.protected<-subset(shark.rich, Protection=="Protected")
shark.rich.unprotected<-subset(shark.rich, Protection=="Unprotected")

#richness of all other species but chondrichthyans
Site.Richness.all<-rowSums(Site.Richness[,which(!colnames(Site.Richness)%in%all.chon)]) # number of species at a site BESIDES chondrichthyans

rich.all.protected<-subset(Site.Richness.all, Protection=="Protected")
rich.all.unprotected<-subset(Site.Richness.all, Protection=="Unprotected")


#identify species found at sites without chondrichthyans

no_shark_site<-Site[Site.SharkAbun==0] # sites that never had any chondrichthyans
Abun_no_sharks<-colSums(subset(Site.Richness, Site%in%no_shark_site)) # number of times a species was found at a site without chondrichthyans
Abun_no_sharks[Abun_no_sharks>0] # species found at sites without chondrichthyans

#number of sites where chondrichthyans have been observed
sum(as.numeric(Site.SharkAbun>0))/167 


#######################################################################
####################Correlations#######################################
#######################################################################

#correlations of chondrichthyan abundance and abundance of other taxa
cor.test(Site.SharkAbun,rowSums(MaxN.all.chon), method="spearman") #correlation with abundance of other taxa
cor(shark.abun.protected, rowSums(MaxN.all.protected), method="spearman") #correlation with abundance of other taxa in protected areas
cor(shark.abun.unprotected, rowSums(MaxN.all.unprotected), method="spearman") #correlation with abundance of other taxa in unprotected areas

#correlations of chondrichthyan abundance and shannon diversity of other taxa
cor.test(Site.SharkAbun,diver.all, method="spearman") #correlation with diversity of other taxa
cor(shark.abun.protected,diver.all.protected) #correlation with diversity of other taxa in protected area
cor(shark.abun.unprotected,diver.all.unprotected) #correlation with diversity of other taxa in unprotected area

#correlations of chondrichthyan richness and abundance of other taxa
cor.test(shark.rich,rowSums(MaxN.all.chon), method="spearman") #correlation with abundance
cor(shark.rich.protected, rowSums(MaxN.all.protected), method="spearman") #correlation with abundance in protected area
cor(shark.rich.unprotected, rowSums(MaxN.all.unprotected), method="spearman") #correlation with abundance in unprotected area

#correlations of chondrichthyan richness and shannon diversity of other taxa
cor.test(shark.rich,diver.all, method="spearman") #correlation with diversity of other taxa
cor(shark.rich.protected,diver.all.protected, method="spearman") #correlation with diversity of taxa in protected area
cor(shark.rich.unprotected,diver.all.unprotected, method="spearman") #correlation with diversity of taxa in unprotected area


#######################################################################
####################GLMs and tests on abundance and diversity##########
#######################################################################
#use gamma distribution since mean Max>0 and continuous
#Do not do pres/abs comparisons for all chondrichthyans due to their high frequency of occurrence
#GLM testing if other taxa more abundant at sites where large sharks observed
l.out1<-glm(rowSums(MaxN.all.large)~factor(large.presab), family=Gamma(link="log")) 
anova(l.out1, update(l.out1, .~-factor(large.presab)), test="Chi") #Likelihood ratio test

#GLM testing if other taxa more abundant at sites where catsharks (shysharks) observed
l.out2<-glm(rowSums(MaxN.all.shyshark)~factor(shyshark.presab),
  family=Gamma(link="log")) 
anova(l.out2, update(l.out2, .~-factor(shyshark.presab)), test="Chi") #Likelihood ratio test

#GLM testing if other taxa more abundant at sites where batoids observed
l.out3<-glm(rowSums(MaxN.all.batoid)~factor(batoid.presab), family=Gamma(link="log")) 
anova(l.out3, update(l.out3, .~-factor(batoid.presab)), test="Chi") #Likelihood ratio test

#t-test effect of large, catshark, and batoid presence on shannon diversity of other taxa
t.test(MaxN.all.shyshark.diver~factor(shyshark.presab), var.equal=TRUE)
t.test(MaxN.all.large.diver~factor(large.presab), var.equal=TRUE)
t.test(MaxN.all.batoid.diver~factor(batoid.presab), var.equal=TRUE)

#check assumption of equal variance:
leveneTest(MaxN.all.shyshark.diver~factor(shyshark.presab))
leveneTest(MaxN.all.large.diver~factor(large.presab))
leveneTest(MaxN.all.batoid.diver~factor(batoid.presab))

#tests if abundance and shannon diversity higher when chondrichthyans are abundance (MaxN > 3)

#GLM testing if other taxa more abundant at sites where batoids observed
l.out_abun<-glm(rowSums(MaxN.site)~factor(manychon.presab), family=Gamma(link="log")) #test if MaxN of all taxa different when chondrichthyans are abundant or not
anova(l.out_abun, update(l.out_abun, .~-factor(manychon.presab)), test="Chi") #Likelihood ratio test

t.test(vegan::diversity(MaxN.site)~factor(manychon.presab), var.equal=TRUE) #test if shannon diversity different when chondrichthyans are abundant or not

leveneTest(vegan::diversity(MaxN.site)~factor(manychon.presab))

#######################################################################
####################PERMANOVAs#########################################
#######################################################################
adonis(MaxN.all.shyshark~shyshark.presab, method="bray", permutations=999) #PERMANOVA for catshark pres/ab effect on other taxa community composition

adonis(MaxN.all.large~large.presab, method="bray", permutations=999) #PERMANOVA for large shark pres/ab effect on other taxa community composition

adonis(MaxN.all.batoid~batoid.presab, method="bray", permutations=999) #PERMANOVA for batoid pres/ab effect on other taxa community composition

adonis(MaxN.all.chon~manychon.presab, method="bray", permutations=999) #PERMANOVA for high chondrichthyan abundance effect on other taxa community composition

################################
#Indicator analysis - DLI's#####
################################
#determine which species are significanty more likely to occur at sites with catsharks (shysharks):
ind.out_shy<-indval(MaxN.all.shyshark, shyshark.presab, numitr=1000) # total DLI for catsharks (species indicating composition when catsharks are present vs absent)
ind.out_shy$maxcls[which(ind.out_shy$pval<=0.05)] #what group/division does the indicator species indicate for (ie large sharks present or absent)?
ind.out_shy$indcls[which(ind.out_shy$pval<=0.05)] #what are the DLI values for each species for the group for which it is the best indicator 

ind.out_large<-indval(MaxN.all.large, large.presab, numitr=1000) # total DLI for species indicating composition when large sharks present vs absent
ind.out_large$maxcls[which(ind.out_large$pval<=0.05)]
ind.out_large$indcls[which(ind.out_large$pval<=0.05)]

ind.out_batoid<-indval(MaxN.all.batoid, batoid.presab, numitr=1000) # total DLI for species indicating composition when batoids present vs absent
ind.out_batoid$maxcls[which(ind.out_batoid$pval<=0.05)]
ind.out_batoid$indcls[which(ind.out_batoid$pval<=0.05)]

ind.out_chon<-indval(MaxN.all.chon, manychon.presab, numitr=1000) # total DLI for species indicating composition when chondrichthyans are abundant vs not abundanct
ind.out_chon$maxcls[which(ind.out_chon$pval<=0.05)]
ind.out_chon$indcls[which(ind.out_chon$pval<=0.05)]

####################################################
#########Reserve analysis###########################
####################################################
#make simulated reserves based on picking a site with or without a mean MaxN of chondrichthyans greater than 3 and the three closest sites to it.
coords_read<-read.csv("SASC_SampleSites.csv", stringsAsFactors=F) #read in locations of our BRUV sites

#mean abundance of species of conservation interest by site
MaxN.interest1<-MaxN[,which(colnames(MaxN)%in%sp.of.interest)]
MaxN.interest<-aggregate(MaxN.interest1, by=list(env$Site, env$Protection), mean)[,-c(1:2)]

coords<-subset(coords_read, coords_read$Station%in%Site) #pick only those sites with BRUVs we dropped in our data set

set.seed(2) #set the seed

#create empty objects to store results of the simulation in a "for loop"
compare.diver.random<-NULL #percent differences for shannon diversity in reserves based on high vs low chondrichthyan abundance
compare.diver.kelp<-NULL #percent differences for shannon diversity in reserves based on high chondrichthyan abundance vs low chondrichthyan abundance on kelp/reef sites

compare.rich.random<-NULL #percent differences for species richness in reserves based on high vs low chondrichthyan abundance
compare.rich.kelp<-NULL #percent differences for species richness in reserves based on high chondrichthyan abundance vs low chondrichthyan abundance on kelp/reef sites

compare.interest.random<-NULL #percent differences for MaxN of species of conservation interest in reserves based on high vs low chondrichthyan abundance
compare.interest.kelp<-NULL #percent differences for MaxN of species of conservation interest in reserves based on high chondrichthyan abundance vs low chondrichthyan abundance on kelp/reef sites

compare.roman.random<-NULL #percent differences for MaxN of red roman in reserves based on high vs low chondrichthyan abundance
compare.roman.kelp<-NULL #percent differences for MaxN of red roman in reserves based on high chondrichthyan abundance vs low chondrichthyan abundance on kelp/reef sites

compare.lobster.random<-NULL #percent differences for MaxN of west coast rock lobster in reserves based on high vs low chondrichthyan abundance
compare.lobster.kelp<-NULL #percent differences for MaxN of west coast rock lobster in reserves based on high chondrichthyan abundance vs low chondrichthyan abundance on kelp/reef sites

mean.roman.random<-NULL #MaxN of red roman in in randomly chosen reserve of low chondrichthyan abundance
mean.roman.kelp<-NULL #MaxN of red roman in in randomly chosen reserve of low chondrichthyan abundance in kelp/reef habitat
mean.roman.shark<-NULL #MaxN of red roman in in randomly chosen reserve of high chondrichthyan abundance

mean.interest.random<-NULL #MaxN of species of conservation interest in in randomly chosen reserve of low chondrichthyan abundance
mean.interest.kelp<-NULL #MaxN of species of conservation interest in in randomly chosen reserve of low chondrichthyan abundance in kelp/reef habitat
mean.interest.shark<-NULL #MaxN of species of conservation interest in in randomly chosen reserve of high chondrichthyan abundance

mean.diver.random<-NULL #shannon diversity in reserve of low chondrichthyan abundance
mean.diver.kelp<-NULL #shannon diversity in reserve of low chondrichthyan abundance in kelp/reef habitat
mean.diver.shark<-NULL #shannon diversity in reserve of high chondrichthyan abundance

mean.rich.random<-NULL #species richness in reserve of low chondrichthyan abundance
mean.rich.kelp<-NULL #species richness in reserve of low chondrichthyan abundance in kelp/reef habitat
mean.rich.shark<-NULL #species richness in reserve of high chondrichthyan abundance

mean.lobster.random<-NULL  #MaxN of west coast rock lobster in in randomly chosen reserve of low chondrichthyan abundance
mean.lobster.kelp<-NULL #MaxN of west coast rock lobster in in randomly chosen reserve of low chondrichthyan abundance in kelp/reef habitat
mean.lobster.shark<-NULL #MaxN of west coast rock lobster in in randomly chosen reserve of high chondrichthyan abundance



for (i in 1:1000) { #repeat simulation 1000 times

chosen_shark<-sample(abundant.shark.sites, 1) #randomly chosen site of high (MaxN>3) chondrichthyan abundance
coords_shark<-subset(coords, Station%in%chosen_shark) #the coordinates of that randomly chosen site of high (MaxN>3) chondrichthyan abundance
n.out_shark<-nn2(coords[,3:4], coords_shark[,3:4], k=4) #3 nearest sites to that randomly chosen site of high (MaxN>3) chondrichthyan abundance
n.out_chosen_shark<-coords[c(n.out_shark$nn.idx),] #name and coordinates of 3 nearest sites to that randomly chosen site of high (MaxN>3) chondrichthyan abundance
mean.diver.shark[i]<-vegan::diversity(colSums(MaxN.site[Site%in%n.out_chosen_shark$Station,])) #shannon diversity of reserve made from 3 nearest sites and that randomly chosen site of high (MaxN>3) chondrichthyan abundance
mean.rich.shark[i]<-sum(as.numeric(colSums(MaxN.site[Site%in%n.out_chosen_shark$Station,])>0)) #species richness of reserve made from 3 nearest sites and that randomly chosen site of high (MaxN>3) chondrichthyan abundance
mean.interest.shark[i]<-sum(as.numeric(colSums(MaxN.interest[Site%in%n.out_chosen_shark$Station,])>0)) #total MaxN (based on summing mean MaxN per site) of species of conservation interest in reserve made from 3 nearest sites and that randomly chosen site of high (MaxN>3) chondrichthyan abundance
mean.roman.shark[i]<-sum(MaxN.interest[Site%in%n.out_chosen_shark$Station,"Red.roman"])  #total MaxN (based on summing mean MaxN per site) of red roman in reserve made from 3 nearest sites and that randomly chosen site of high (MaxN>3) chondrichthyan abundance
mean.lobster.shark[i]<-sum(MaxN.interest[Site%in%n.out_chosen_shark$Station,"West.coast.rock.lobster"]) #total MaxN (based on summing mean MaxN per site) of west coast rock lobster in reserve made from 3 nearest sites and that randomly chosen site of high (MaxN>3) chondrichthyan abundance

chosen_random<-sample(Sites_low_sharks, 1) #randomly chosen site of low (MaxN<=3) chondrichthyan abundance
coords_random<-subset(coords, Station%in%chosen_random)  #the coordinates of that randomly chosen site of low (MaxN<=3) chondrichthyan abundance
n.out_random<-nn2(coords[,3:4], coords_random[,3:4], k=4) #3 nearest sites to that randomly chosen site of low (MaxN<=3) chondrichthyan abundance
n.out_chosen_random<-coords[c(n.out_random$nn.idx),] #name and coordinates of 3 nearest sites to that randomly chosen site of low (MaxN<=3) chondrichthyan abundance
mean.diver.random[i]<-vegan::diversity(colSums(MaxN.site[Site%in%n.out_chosen_random$Station,])) #shannon diversity of reserve made from these four sites
compare.diver.random[i]<-(mean.diver.shark[i]-mean.diver.random[i])/mean(c(mean.diver.shark[i],mean.diver.random[i]))  #compare % difference in shannon diversity of reserves made based on high vs low chondrichthyan abundance
mean.rich.random[i]<-sum(as.numeric(colSums(MaxN.site[Site%in%n.out_chosen_random$Station,])>0)) #total species richness of reserve made from these four sites
compare.rich.random[i]<-(mean.rich.shark[i]-mean.rich.random[i])/mean(c(mean.rich.shark[i],mean.rich.random[i])) #compare % difference in species richness of reserves made based on high vs low chondrichthyan abundance
mean.interest.random[i]<-sum(as.numeric(colSums(MaxN.interest[Site%in%n.out_chosen_random$Station,])>0)) #total MaxN (based on summing mean MaxN per site) for species of conservation interest of reserve made from these four sites
compare.interest.random[i]<-(mean.interest.shark[i]-mean.interest.random[i])/mean(c(mean.interest.shark[i],mean.interest.random[i])) #compare % difference in MaxN of species of conservation interest of reserves made based on high vs low chondrichthyan abundance
mean.roman.random[i]<-sum(MaxN.interest[Site%in%n.out_chosen_random$Station,"Red.roman"]) #total MaxN (based on summing mean MaxN per site) for red roman of reserve made from these four sites
compare.roman.random[i]<-(mean.roman.shark[i]-mean.roman.random[i])/mean(c(mean.roman.shark[i],mean.roman.random[i])) #compare % difference in MaxN of red roman of reserves made based on high vs low chondrichthyan abundance
mean.lobster.random[i]<-sum(MaxN.interest[Site%in%n.out_chosen_random$Station,"West.coast.rock.lobster"]) #total MaxN (based on summing mean MaxN per site) for west coast rock lobster of reserve made from these four sites
compare.lobster.random[i]<-(mean.lobster.shark[i]-mean.lobster.random[i])/mean(c(mean.lobster.shark[i],mean.lobster.random[i])) #compare % difference in MaxN of west coast rock lobster of reserves made based on high vs low chondrichthyan abundance

chosen_kelp<-sample(ReefKelp_low_sharks_sites, 1) # random kelp or reef site with low (MaxN<=3) mean chondrichthyan abundance
coords_kelp<-subset(coords, Station%in%chosen_kelp)  #the coordinates of that randomly chosen kelp/reef site of low (MaxN<=3) chondrichthyan abundance
n.out_kelp<-nn2(coords[,3:4], coords_kelp[,3:4], k=4) #3 nearest sites to that randomly chosen kelp/reef site of low (MaxN<=3) chondrichthyan abundance
n.out_chosen_kelp<-coords[c(n.out_kelp$nn.idx),] #name and coordinates of 3 nearest sites to that randomly chosen kelp/reef site of low (MaxN<=3) chondrichthyan abundance
mean.diver.kelp[i]<-vegan::diversity(colSums(MaxN.site[Site%in%n.out_chosen_kelp$Station,])) #shannon diversity of reserve made from these four sites
compare.diver.kelp[i]<-(mean.diver.shark[i]-mean.diver.kelp[i])/mean(c(mean.diver.shark[i],mean.diver.kelp[i]))  #compare % difference in shannon diversity of reserves made based on high vs low chondrichthyan abundance
mean.rich.kelp[i]<-sum(as.numeric(colSums(MaxN.site[Site%in%n.out_chosen_kelp$Station,])>0))  #total species richness of reserve made from these four sites
compare.rich.kelp[i]<-(mean.rich.shark[i]-mean.rich.kelp[i])/mean(c(mean.rich.shark[i],mean.rich.kelp[i])) #compare % difference in species richness of reserves made based on high vs low chondrichthyan abundance
mean.interest.kelp[i]<-sum(as.numeric(colSums(MaxN.interest[Site%in%n.out_chosen_kelp$Station,])>0)) #total MaxN (based on summing mean MaxN per site) for species of conservation interest of reserve made from these four sites
compare.interest.kelp[i]<-(mean.interest.shark[i]-mean.interest.kelp[i])/mean(c(mean.interest.shark[i],mean.interest.kelp[i]))  #compare % difference in MaxN of species of conservation interest of reserves made based on high vs low chondrichthyan abundance
mean.roman.kelp[i]<-sum(MaxN.interest[Site%in%n.out_chosen_kelp$Station,"Red.roman"]) #total MaxN (based on summing mean MaxN per site) for red roman of reserve made from these four sites
compare.roman.kelp[i]<-(mean.roman.shark[i]-mean.roman.kelp[i])/mean(c(mean.roman.shark[i],mean.roman.kelp[i]))  #compare % difference in MaxN of red roman of reserves made based on high vs low chondrichthyan abundance
mean.lobster.kelp[i]<-sum(MaxN.interest[Site%in%n.out_chosen_kelp$Station,"West.coast.rock.lobster"]) #total MaxN (based on summing mean MaxN per site) for west coast rock lobster of reserve made from these four sites
compare.lobster.kelp[i]<-(mean.lobster.shark[i]-mean.lobster.kelp[i])/mean(c(mean.lobster.shark[i],mean.lobster.kelp[i])) #compare % difference in MaxN of west coast rock lobster of reserves made based on high vs low chondrichthyan abundance


}

#mean % difference between reserves based on sites with lots or few chondrichthyans (including those restricted to be in kelp/reef habitat of low chondrichthyan abundance)

#for MaxN of red romans:
mean(compare.roman.random*100)
mean(compare.roman.kelp*100)

qnorm(0.975)*sd(compare.roman.random*100)/sqrt(length(compare.roman.random)) #95% CI
qnorm(0.975)*sd(compare.roman.kelp*100)/sqrt(length(compare.roman.kelp))  #95% CI

#comparing for MaxN of west coast rock lobster
compare.lobster.random<-compare.lobster.random[!is.na(compare.lobster.random)] #remove NAs (no rock lobsters in either site)
compare.lobster.kelp<-compare.lobster.kelp[!is.na(compare.lobster.kelp)] #remove NAs (no rock lobsters in either site)

mean(compare.lobster.random*100)
mean(compare.lobster.kelp*100)

qnorm(0.975)*sd(compare.lobster.random*100)/sqrt(length(compare.lobster.random))  #95% CI
qnorm(0.975)*sd(compare.lobster.kelp*100)/sqrt(length(compare.lobster.kelp))  #95% CI

#comparing MaxN for species of conservation concern
mean(compare.interest.random*100)
mean(compare.interest.kelp*100)

qnorm(0.975)*sd(compare.interest.random*100)/sqrt(length(compare.interest.random)) #95% CI
qnorm(0.975)*sd(compare.interest.kelp*100)/sqrt(length(compare.interest.kelp)) #95% CI

#comparing overall shannon diversity
mean(compare.diver.random*100)
mean(compare.diver.kelp*100)
qnorm(0.975)*sd(compare.diver.random*100)/sqrt(length(compare.diver.random)) #95% CI
qnorm(0.975)*sd(compare.diver.kelp*100)/sqrt(length(compare.diver.kelp)) #95% CI

#compare species richness
mean(compare.rich.random*100)
mean(compare.rich.kelp*100)
qnorm(0.975)*sd(compare.rich.random*100)/sqrt(length(compare.rich.random)) #95% CI
qnorm(0.975)*sd(compare.rich.kelp*100)/sqrt(length(compare.rich.kelp)) #95% CI

  