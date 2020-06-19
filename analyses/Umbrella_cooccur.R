##########################
##Co-occurrence analysis##
##########################

library(cooccur) # co-occurrence analysis
library(igraph) #centrality analysis
library(CINNA) #centrality analysis
library(sna) #centrality analysis
library(dplyr) #for %>% function and distinct function
library(tibble) # for rowid_to_column function

#Reading in data
BRUV.MaxN<-read.csv("data/Osgood_etal_2020_MaxN_Data.csv", header=TRUE, stringsAsFactors=FALSE)
env<-read.csv("data/Osgood_etal_2020_BRUV_sampling_metadata_environmentals.csv", header=TRUE, stringsAsFactors=FALSE)

########################################################################
##Remove species (ie invertebrates) that are not relevant/ are rare#####
########################################################################
##Remove species (ie invertebrates) that are not relevant/ are rare
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


#Setting the names of various ways of breaking up chondrichthyan class so can subset later - need periods as R deletes the spaces
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
          "Spotted.gully.shark",
               "Bronze.whaler", "Smooth.hammerhead.shark", 
               "Smooth.hound.shark",  
                "Soupfin.shark")

all.batoid<-c("Lesser.guitarfish", "Eagle.ray", "Biscuit.skate", 
  "Spearnose.skate", "Short.tail.stingray")

all.crust<-c("Cape.rock.crab", "Hermit.crab", "Masked.crab", "Masked.crab.sp2", 
  "Sandflat.crab", "Three.spotted.swimming.crab", "West.coast.rock.lobster")

all.ceph<-c("Common.octopus", "Chokka.squid")

all.teleost<-c("Bank.steenbras", "Barehead.goby","Barred.fingerfin", "Beaked.sandfish",
  "Black.seacatfish", "Blacktail", "Blue.hottentot", "Bluefin.gurnard", "Cape.horse.mackerel",
      "Cape.sole", "Carpenter", "Evileye.blaasop", "Fransmadam", "Galjoen", "Geelbek", 
      "Giant.yellowtail", "Hottentot", "Janbruin", "Jutjaw", "Longsnout.pipefish", "Panga",
      "Red.roman", "Red.steenbras", "Red.stumpnose", "Redfingers", "Shad", "Silver.kob",
      "Slender.baardman", "Snakelet", "Spinynose.horsefish", "Steentjie", "Strepie", 
      "Super.klipfish", "Twotone.fingerfin", "White.seacatfish", "White.stumpnose", 
      "Yellowback.fusilier", "Yellowbelly.rock.cod", "Zebra")

sp.of.interest<-c("Red.steenbras", "Galjoen", "Geelbek", "Red.roman", 
  "Red.stumpnose", "Janbruin", "Blacktail", "Shad", "Silver.kob", "Slender.baardman",
  "West.coast.rock.lobster", "White.stumpnose", "Yellowbelly.rock.cod", "Zebra")

MaxN<-BRUV.MaxN[,-c(1:5)] # Creates abundance only data frame (remove "meta-data" columns)
#head(MaxN) #Check that I was successful by looking at first few rows

BRUV.Richness<-apply(MaxN, 2, function(x) {as.numeric(x>0)}) # Create data frame that is presence/absence (0 and 1) only


#Use below for co-occurrence by site rather than BRUV drop:

Site.Richness<-aggregate(BRUV.Richness, by=list(env$Site), sum)[,-1] # total number of times each species found per site
Site.Richness<-apply(Site.Richness, 2, function(x) {as.numeric(x>0)}) # Returns data frame to 0 and 1 only - ie yes if species EVER found there, no if NEVER found there

#######################################
#####Co-occurrence analysis############
#######################################

#co-occurrrence values for all species:
co.out.Site<-cooccur(t(Site.Richness), spp_names=TRUE, type="spp_site", thresh=FALSE) # Do co-occur analysis based on site level


#results of co-occurrence analysis into data frame that can  be worked with:
print.site<-co.out.Site$results


#manually calculate effect sizes (difference between observed and expected co-occurrence standardized by number of sites):
print.site$effects<-(print.site$obs_cooccur-print.site$exp_cooccur)/167

print.site_pos<-subset(print.site, p_gt<0.05) # significant positive co-occurrences

#use this to assess the 20 strongest interactions (first twenty rows):
print.site_pos<-print.site_pos[order(-print.site_pos$effects),] # significant positive co-occurrences ordered by effect size
print.site_pos[1:20,]
#subsetting co-occurrence results for interactions involving chondrichthyans:
all.results_chon<-subset(print.site, sp1_name%in%all.chon | sp2_name%in%all.chon) 

sig_pos_chon<-subset(all.results_chon, p_gt<0.05) # significant positive co-occurrences with chondrichthyans
sig_neg_chon<-subset(all.results_chon, p_lt<0.05) # significant negative co-occurrences with chondrichthyans
nrow(sig_pos_chon)/nrow(all.results_chon) #percent that are significant positive associations with chondrichthyans
nrow(sig_neg_chon)/nrow(all.results_chon) #percent that are significant negative associations with chondrichthyans

#subsetting co-occurrence results for interactions involving chondrichthyans to species of conservation interest:

all.results_int1<-subset(print.site, sp1_name%in%all.chon & sp2_name%in%sp.of.interest) #chondrichthyans are in the first column
all.results_int2<-subset(print.site, sp1_name%in%sp.of.interest & sp2_name%in%all.chon) #chondrichthyans are in the second column

all.results_int<-rbind(all.results_int1, all.results_int2) # get data frame of all interactions of chondrichthyans to species of conservation interest (regardless if they were in first or second column)

sig_pos_int<-subset(all.results_int, p_gt<0.05) #subset to just significant positive associations

sig_neg_int<-subset(all.results_int, p_lt<0.05) #subset to just significant negative associations

nrow(sig_pos_int)/nrow(all.results_int) #percent of significant positive associations
nrow(sig_neg_int)/nrow(all.results_int) #percent of significant negative associations


#subsetting co-occurrence results for all interactions involving species of conservation interest:
all.results_interest<-subset(print.site, sp1_name%in%sp.of.interest | sp2_name%in%sp.of.interest)

sig_pos_interest<-subset(all.results_interest, p_gt<0.05) #subset to just significant positive associations
sig_neg_interest<-subset(all.results_interest, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_interest) #number of significant positive associations
nrow(sig_neg_interest) #number of significant negative associations


#subsetting co-occurrence results for all interactions involving teleosts:
all.results_tel<-subset(print.site, sp1_name%in%all.teleost | sp2_name%in%all.teleost)

sig_pos_tel<-subset(all.results_tel, p_gt<0.05) #subset to just significant positive associations
sig_neg_tel<-subset(all.results_tel, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_tel)/nrow(all.results_tel) #number of significant positive associations
nrow(sig_neg_tel)/nrow(all.results_tel) #number of significant negative associations


#subsetting co-occurrence results for all interactions involving cephalopods:
all.results_ceph<-subset(print.site, sp1_name%in%all.ceph | sp2_name%in%all.ceph)

sig_pos_ceph<-subset(all.results_ceph, p_gt<0.05) #subset to just significant positive associations
sig_neg_ceph<-subset(all.results_ceph, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_ceph) #number of significant positive associations
nrow(sig_neg_ceph) #number of significant negative associations


#subsetting co-occurrence results for all interactions involving crustaceans:
all.results_crust<-subset(print.site, sp1_name%in%all.crust | sp2_name%in%all.crust)

sig_pos_crust<-subset(all.results_crust, p_gt<0.05) #subset to just significant positive associations
sig_neg_crust<-subset(all.results_crust, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_crust) #number of significant positive associations
nrow(sig_neg_crust) #number of significant negative associations


#subsetting co-occurrence results for all interactions involving catsharks (same as shysharks):
all.results_shy<-subset(print.site, sp1_name%in%all.shyshark | sp2_name%in%all.shyshark)

sig_pos_shy<-subset(all.results_shy, p_gt<0.05) #subset to just significant positive associations
sig_neg_shy<-subset(all.results_shy, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_shy)/nrow(all.results_shy) #percent that are significant positive associations
nrow(sig_neg_shy)/nrow(all.results_shy) #percent that are significant negative associations


#subsetting co-occurrence results for all interactions involving large sharks:
all.results_large<-subset(print.site, sp1_name%in%all.large.shark | sp2_name%in%all.large.shark)

sig_pos_large<-subset(all.results_large, p_gt<0.05) #subset to just significant positive associations
sig_neg_large<-subset(all.results_large, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_large) #number of significant positive associations
nrow(sig_neg_large) #number of significant negative associations

#individual associations of sevengill shark:
sevengill_associations<-subset(all.results_large, sp1_name=="Broadnose.sevengill.shark" | sp2_name=="Broadnose.sevengill.shark")

#subsetting co-occurrence results for all interactions involving large sharks:
all.results_batoid<-subset(print.site, sp1_name%in%all.batoid| sp2_name%in%all.batoid)

sig_pos_batoid<-subset(all.results_batoid, p_gt<0.05) #subset to just significant positive associations
sig_neg_batoid<-subset(all.results_batoid, p_lt<0.05) #subset to just significant negative associations
nrow(sig_pos_batoid) #number of significant positive associations
nrow(sig_neg_batoid) #number of significant negative associations

#species with co-occurring strongly with roman and west coast rock lobster:
roman_lobster_results<-subset(print.site, sp1_name%in%c("Red.roman", "West.coast.rock.lobster") | sp2_name%in%c("Red.roman", "West.coast.rock.lobster"))
roman_lobster_results<-roman_lobster_results[order(-roman_lobster_results$effects),]
head(roman_lobster_results)

######
##Examining how co-occurrence to this species of interest changes based on protection:
#break BRUV.Richness by habitat

protected.richness<-BRUV.Richness[which(env$Protection=="Protected"),] #sites in protected areas
unprotected.richness<-BRUV.Richness[which(env$Protection=="Unprotected"),] #sites outside of protected areas

env.protected<-env[which(env$Protection=="Protected"),] #get list of sites in protected areas
env.unprotected<-env[which(env$Protection=="Unprotected"),] #get list of sites outside of protected areas

#aggregated by site inside protected areas:
Site.Richness.Protected1<-aggregate(protected.richness, by=list(env.protected$Site), sum)[,-1] # total number of times each species found per site
Site.Richness.Protected<-apply(Site.Richness.Protected1, 2, function(x) {as.numeric(x>0)}) # Returns data frame to 0 and 1 only - ie yes if species EVER found there, no if NEVER found there

#aggregated by site outside protected areas:
Site.Richness.Unprotected1<-aggregate(unprotected.richness, by=list(env.unprotected$Site), sum)[,-1] # total number of times each species found per site
Site.Richness.Unprotected<-apply(Site.Richness.Unprotected1, 2, function(x) {as.numeric(x>0)}) # Returns data frame to 0 and 1 only - ie yes if species EVER found there, no if NEVER found there

#co-occurrence inside protected areas:
co.out.protected<-cooccur(t(Site.Richness.Protected), spp_names=TRUE, type="spp_site", thresh=TRUE) # Do co-occur analysis based on site level
print.site.protected<-co.out.protected$results # get list of interactions into its own object so I can subset to focus on roman and lobster and catsharks
print.site.protected$effects<-(print.site.protected$obs_cooccur-print.site.protected$exp_cooccur)/nrow(Site.Richness.Protected) #calculate effect sizes (dividing by number of protected sites)
all.results_protected<-subset(print.site.protected, sp1_name%in%all.chon | sp2_name%in%all.chon) #subset to just chondrichthyans
all.results_protected_roman_lobster<-subset(all.results_protected, sp1_name%in%c("Red.roman", "West.coast.rock.lobster") | sp2_name%in%c("Red.roman", "West.coast.rock.lobster")) #subset to just red roman and rock lobster


#co-occurrence outside protected areas:
co.out.unprotected<-cooccur(t(Site.Richness.Unprotected), spp_names=TRUE, type="spp_site", thresh=TRUE) # Do co-occur analysis based on site level
print.site.unprotected<-co.out.unprotected$results # get list of interactions into its own object so I can subset to focus on roman and lobster and catsharks
print.site.unprotected$effects<-(print.site.unprotected$obs_cooccur-print.site.unprotected$exp_cooccur)/nrow(Site.Richness.Unprotected) #calculate effect sizes (dividing by number of protected sites)
all.results_unprotected<-subset(print.site.unprotected, sp1_name%in%all.chon | sp2_name%in%all.chon) #subset to just chondrichthyans
all.results_unprotected_roman_lobster<-subset(all.results_unprotected, sp1_name%in%c("Red.roman", "West.coast.rock.lobster") | sp2_name%in%c("Red.roman", "West.coast.rock.lobster")) #subset to just red roman and rock lobster

#######

####################################################################################
#####################Co-occurrences by species for each group#######################
####################################################################################

par.out_site<-pair.attributes(co.out.Site)# get summary of co-occurrences for each species by species

####################################################
####for chondrichthyans#############################
####################################################
par.out.chon_site<-par.out_site[which(par.out_site$sppname%in%all.chon),] #Summary of co-occurrences for chondrichthyan species 

par.out_site.shy<-par.out_site[which(par.out_site$sppname%in%all.shyshark),] #Selecting catsharks for example

par.out_site.batoid<-par.out_site[which(par.out_site$sppname%in%all.batoid),] #Selecting all batoids

par.out_site.large<-par.out_site[which(par.out_site$sppname%in%all.large.shark),] #Selecting all large sharks

#for teleosts
par.out.teleost_site<-par.out_site[which(par.out_site$sppname%in%all.teleost),] #Selecting all teleosts

#species of interest
par.out.interest_site<-par.out_site[which(par.out_site$sppname%in%sp.of.interest),] #Selecting species of conservation interest

#crustaceans
par.out.crust_site<-par.out_site[which(par.out_site$sppname%in%all.crust),] #Selecting all crustaceans

#cephalopods
par.out.ceph_site<-par.out_site[which(par.out_site$sppname%in%all.ceph),] #Selecting all cephalopods

##########################################################################
###Mean percent positive and negative co-occurrences by species groups: ##
##########################################################################
#Below is the mean percent of positive, negative, and random co-occurrences for a species group (the mean % is averaged across all species in that group)

#Chondrichthyans
mean_chondr_pos<-mean(par.out.chon_site$pos)
mean_chondr_neg<-mean(par.out.chon_site$neg)
mean_chondr_random<-mean(par.out.chon_site$rand)

#catsharks
mean_shy_pos<-mean(par.out_site.shy$pos)
mean_shy_neg<-mean(par.out_site.shy$neg)
mean_shy_random<-mean(par.out_site.shy$rand)

#large sharks
mean_large_pos<-mean(par.out_site.large$pos)
mean_large_neg<-mean(par.out_site.large$neg)
mean_large_random<-mean(par.out_site.large$rand)

#batoids
mean_batoid_pos<-mean(par.out_site.batoid$pos)
mean_batoid_neg<-mean(par.out_site.batoid$neg)
mean_batoid_random<-mean(par.out_site.batoid$rand)

#teleosts
mean_tel_pos<-mean(par.out.teleost_site$pos)
mean_tel_neg<-mean(par.out.teleost_site$neg)
mean_tel_random<-mean(par.out.teleost_site$rand)

#crustaceans
mean_crust_pos<-mean(par.out.crust_site$pos)
mean_crust_neg<-mean(par.out.crust_site$neg)
mean_crust_random<-mean(par.out.crust_site$rand)

#cephalopods
mean_ceph_pos<-mean(par.out.ceph_site$pos)
mean_ceph_neg<-mean(par.out.ceph_site$neg)
mean_ceph_random<-mean(par.out.ceph_site$rand)

####################################
##Comparing mean effect sizes: #####
####################################

#chondrichthyan effect size means and sd positive:
mean(sig_pos_chon$effects)
sd(sig_pos_chon$effects)
mean(sig_pos_interest$effects)
mean(sig_pos_tel$effects)

#chondrichthyan effect size means and sd negative:
mean(sig_neg_chon$effects)
sd(sig_neg_chon$effects)
mean(sig_neg_interest$effects)
mean(sig_neg_tel$effects)

##mean effect sizes per chondrichthyan species:

#data frame of positive effect sizes for each individual species interaction:
Chon_effects_pos<-data.frame(species=c(as.character(sig_pos_chon$sp1_name), 
  as.character(sig_pos_chon$sp2_name)), 
  effects=c(sig_pos_chon$effects, sig_pos_chon$effects))
Chon_effects_pos<-subset(Chon_effects_pos, species%in%all.chon)

#data frame of negative effect sizes for each individual species interaction:
Chon_effects_neg<-data.frame(species=c(as.character(sig_neg_chon$sp1_name), 
  as.character(sig_neg_chon$sp2_name)), 
  effects=c(sig_neg_chon$effects, sig_neg_chon$effects))
Chon_effects_neg<-subset(Chon_effects_neg, species%in%all.chon)

pos_effect_Chon<-aggregate(effects~species, Chon_effects_pos, mean) # mean positive effect size per species across all their interactions
neg_effect_Chon<-aggregate(effects~species, Chon_effects_neg, mean) # mean negative effect size per species across all their interactions


##mean effect sizes per teleost species:

#data frame of positive effect sizes for each individual teleost species interaction:
Tel_effects_pos<-data.frame(species=c(as.character(sig_pos_tel$sp1_name), 
  as.character(sig_pos_tel$sp2_name)), 
  effects=c(sig_pos_tel$effects, sig_pos_tel$effects))
Tel_effects_pos<-subset(Tel_effects_pos, species%in%all.teleost)

#data frame of negative effect sizes for each individual species teleost interaction:
Tel_effects_neg<-data.frame(species=c(as.character(sig_neg_tel$sp1_name), 
  as.character(sig_neg_tel$sp2_name)), 
  effects=c(sig_neg_tel$effects, sig_neg_tel$effects))
Tel_effects_neg<-subset(Tel_effects_neg, species%in%all.teleost)

pos_effect_tel<-aggregate(effects~species, Tel_effects_pos, mean) # mean positive effect size per species across all their interactions
neg_effect_tel<-aggregate(effects~species, Tel_effects_neg, mean) # mean negative effect size per species across all their interactions

################################################
####Friendship paradox##########################
################################################

#determine for how many species each species was one of the top three most strongly co-occurring positively
#subset by each species and save the identity of the top three species that co-occur most strongly and positively with it
#loop over all the species to identify the top three most strongly co-ocurring species for each species

MaxN_subset<-MaxN[,which(colSums(MaxN)>1)] #exclude species who occur only once

trial_names<-colnames(MaxN_subset) # save the names for use in subsetting and for loop (loop over the species names)

#set up the initial entries (first three) of each vector to store the results:
temp_cooccur<-subset(print.site, sp1_name%in%trial_names[1] | sp2_name%in%trial_names[1]) #first species to assess which species co-occured most strongly with it (top 3)
 #make the name character for when the species is both in column 1 (sp1) and column 2 (sp2) of the print.site co-occurr results
temp_cooccur$sp1_name<-as.character(temp_cooccur$sp1_name)
temp_cooccur$sp2_name<-as.character(temp_cooccur$sp2_name)

temp_cooccur<-temp_cooccur[order(-(temp_cooccur$obs_cooccur-temp_cooccur$exp_cooccur)),] #order interactions for that species by strongest interactions/effects
result_temp<-c(temp_cooccur$sp1_name[1:3], temp_cooccur$sp2_name[1:3]) #identify the three strongest species that interact with it
result<-result_temp[result_temp!=trial_names[1]] #get the names of the species that co-occur with the target species (not the name of the target species itself)
#the code returns species in both column 1 and column 2, but one of those will be the target species. I just want the names of the co-occurring species
#Since African penguin is first alphabetically, the first three entries represent the three top species most strongly interacting with the penguin.

#do the same thing as above but now for the other species:
#run the for loop on the remaining rows
for (i in 2:length(trial_names)) {


temp_cooccur<-subset(print.site, sp1_name%in%trial_names[i] | sp2_name%in%trial_names[i])
temp_cooccur$sp1_name<-as.character(temp_cooccur$sp1_name)
temp_cooccur$sp2_name<-as.character(temp_cooccur$sp2_name)

temp_cooccur<-temp_cooccur[order(-(temp_cooccur$obs_cooccur-temp_cooccur$exp_cooccur)),]
result_temp<-c(temp_cooccur$sp1_name[1:3], temp_cooccur$sp2_name[1:3])
result<-c(result, result_temp[result_temp!=trial_names[i]])

}

t.out<-table(result) #see how many times each species was one of the top three most strongly interacting species
t.out[order(-t.out)] #order the table from most to least times as one of the top three


####################################
##Network and centrality analysis###
####################################

#print.site is my list of co-occuring species

co_occurrences<-subset(print.site, obs_cooccur>0) #only for species that did co-occurr with other things

#get list of distinct species from both columns of co-occur results:
#column 1:
sp1 <- co_occurrences %>%
  distinct(sp1_name) %>%
  rename(label = sp1_name)

#column 2:
sp2 <- co_occurrences %>%
  distinct(sp2_name) %>%
  rename(label = sp2_name)

nodes <- full_join(sp1, sp2, by = "label") %>% rowid_to_column("id") #combine my two lists of uniques species from each column
nodes$id<-as.character(nodes$id) #make the id column character (going to use numerical id to make the "graphs"/networks)
nodes$label<-as.character(nodes$label) #ensure the name labels are also characters rather than factors

#ensure the co-occurrence data frame is ordered by species in first column then second
per_network <- co_occurrences %>%  
  group_by(sp1_name, sp2_name) %>%
  ungroup()

#addes "from" and "to" columns necessary to make the network (ie interaction from this species to this other species - from and to just indicated the network, direction will not matter)
edges <- per_network %>% 
  left_join(nodes, by = c("sp1_name" = "label")) %>% 
  rename(from = id) %>% 
  left_join(nodes, by = c("sp2_name" = "label")) %>% 
  rename(to = id)

#select just from and to columns (the essentials of the network)
edges <- select(edges, from, to)

#make a "graph" (ie network) from the links defined by from and to above:
shark_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

#Centrality measures from using igraph:

centrality.out<-cbind(edges, edge_betweenness(shark_igraph)) #Betweenness centrality
names(centrality.out)<-c("from", "to", "centrality") #better names
centrality.out<-centrality.out[order(-centrality.out$centrality),] #order by betweeness centrality
head(centrality.out) #see the first six most central by this measure

eigen.out<-cbind(nodes, eigen_centrality(shark_igraph)$vector) #eigenvalue centrality
names(eigen.out)<-c("id", "label", "centrality") #better names
eigen.out<-eigen.out[order(-eigen.out$centrality),] #order by eigenvalue centrality
head(eigen.out) #see the first six most central by this measure

degree.out<-cbind(nodes, graph.strength(shark_igraph)) #degree centrality
names(degree.out)<-c("id", "label", "centrality") #better names
degree.out<-degree.out[order(-degree.out$centrality),] #order by degree centrality
head(degree.out) #see the first six most central by this measure


#A faster way to summarize all the centrality measures in a matrix of species by centrality measure (from three packages:
  res <- matrix(0,vcount(shark_igraph),12) #create matrix (12 columns for 12 centrality measures)
 
 #degree, betweenness, closeness, and eigenvalue centrality calculated using igraph package:
  res[,1] <- igraph::degree(shark_igraph)
  res[,2] <- igraph::betweenness(shark_igraph)
  res[,3] <- igraph::closeness(shark_igraph)
  res[,4] <- igraph::eigen_centrality(shark_igraph)$vector
  
  A <- get.adjacency(shark_igraph,sparse=F) #modified graph into class for use in sna package

 #five other centrality measures from the sna package:
  res[,5] <- sna::flowbet(A) 
  res[,6] <- sna::loadcent(A) 
  res[,7] <- sna::gilschmidt(A) 
  res[,8] <- sna::infocent(A) 
  res[,9] <- sna::stresscent(A) 

 #three other centrality measures from the CINNA package:
  res[,10] <- CINNA::dangalchev_closeness_centrality(shark_igraph)
  res[,11] <- CINNA::harmonic_centrality(shark_igraph)
  res[,12] <- 1/CINNA::local_bridging_centrality(shark_igraph) 


 res1 <-  apply(res,2,function(x) round(x,8)) #round the centrality measures

which.max<- function(x) which(max(x)==x)

apply(res1, 2, which.max)

#see the top six species for select centrality measures, once that measure has been ordered from most to least central species:
head(cbind(nodes,res1[,1])[order(-res1[,1]),])
head(cbind(nodes,res1[,2])[order(-res1[,2]),])
head(cbind(nodes,res1[,12])[order(-res1[,12]),])

head(cbind(nodes,res1[,3])[order(-res1[,3]),])
head(cbind(nodes,res1[,4])[order(-res1[,4]),])




