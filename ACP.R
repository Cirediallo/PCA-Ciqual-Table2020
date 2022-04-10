install.packages("VIM")
install.packages("funModeling") 
install.packages("summarytools")
library(funModeling)
library(VIM)

Data = read.csv("C:/Users/Mamadou/Documents/Cours/M1 ATAL/Analyse de données/projet/Analysis/Table Ciqual 2020_FR_2020 07 07SsInf0.csv", sep=";")
Data
str(Data)
API = subset(Data, Data$alim_ssgrp_nom_fr=="fromages et assimilés")
API
dim(API)
typeof(API)
is.data.frame(API)
#Df = API  %>% select(-(alim_grp_code:alim_nom_sci))
Frame = subset(API, select=(-c(alim_grp_code:alim_nom_sci)))

for ( iter in 1:ncol(Frame)){Frame[,iter] = as.numeric(gsub(",", ".", gsub("\\.", "", Frame[,iter])))}
FrameWithName = cbind(API$alim_nom_fr, API$alim_ssssgrp_nom_fr, API$alim_code, Frame) #we add useful columns
str(FrameWithName) #liste les variables, indique leur type ainsi qu'un echantillon des 1ères valeurs
NA_Studies = aggr(FrameWithName,
                  col=c('navyblue','red'),
                  numbers=TRUE,
                  sortVars=TRUE,
                  labels=names(data),
                  cex.axis=.4, gap=3,
                  ylab=c("Histogram of missing data","Pattern"))

summary(FrameWithName)

df_status(Frame) #connaitre le nombre de données manquantes et le pourcentage correspondant pour chaque variable

#columnToDelete = subset(Frame, df_status(Frame)$p_na > 10)
#remove column with more than 10 NA's values which represent 
removenacol = function(mat){
  nacols = colSums(is.na(mat))
  for(i in ncol(mat):1){
    if(nacols[i] > 10){
      mat = mat[-c(i)]
    }
  }
  return(mat)
}




newDataFrame = removenacol(Frame)

newDataFrame
dim(newDataFrame)
#====================================================================================#
#==== A CE NIVEAU DECIDER DE SUPPRIMER TOUTES LES VALEURS SUPERIEURES A 8% ========
#====================================================================================#

#Cette fonction est dans le rapport à citer mais à ne pas mettre dans le code
summarytools::descr(FrameWithName,
                    headings = FALSE,
                    transpose = TRUE)

#Analyse UNIDIMENSIONNELS
#Graphique: permettent de savoir si les suppositions(faites sur les variables) sont plus ou moins bien respectés c'est ce qu'on appelle analyse exploratoire des données EDA(Exploratory Data ANalysis) en anglais

boxplot(Frame$Energie..Règlement.UE.N..1169.2011..kJ.100.g., main="Analyse de la distribution de la variable Energie Règlement UE N1169")
boxplot(Frame$Beta.Carotène..µg.100.g.) #dissimetrie avec un <<fort>> étalement vers les grandes valeurs
boxplot(Frame$Sélénium..µg.100.g.)
boxplot(Frame$AG.polyinsaturés..g.100.g., main="Distribution de la variable AG Polyinsaturés")
boxplot(Frame$AG.monoinsaturés..g.100.g., main="Analyse de la distribution de la variable AG Monoinsaturés")
boxplot(Frame)




cleanData = na.omit(newDataFrame)
View(cleanData)
dim(cleanData) #maintenant on travaille avec une matrice de 106 lignes et 15 variables


APICleaned = cleanData[-1,] #suppression de la 1ère ligne qui correspond aux fromages moyens
dim(APICleaned)
View(APICleaned)


# ==================================================================== #
# ************************ PART 2 *********************************** #
# ********************** ANALYSE UNIDIMENSIONNELLE ****************** #
# =================================================================== #

summary(APICleaned)
describe(APICleaned) # a ce niveau on voit avec la colonne distinct  qu'on a des doublon et on peut voir la fréquence d'apparition avec freq
freq(APICleaned$Eau..g.100.g.)
plot(density(APICleaned$Eau..g.100.g.), main = "Eau")

boxplot(APICleaned, main="Ensemble des données", ylab="valeurs")
boxplot(APICleaned$Glucides..g.100.g.)
boxplot(APICleaned$AG.monoinsaturés..g.100.g.) #distribution proche de la normale
boxplot(Frame$AG.monoinsaturés..g.100.g., main="Distribution de la variable AG Monoinsaturés")
boxplot(APICleaned$AG.polyinsaturés..g.100.g., main="Distribution polyinsaturés", ylab="valeurs")
boxplot(APICleaned$Eau..g.100.g., main= "Distribution Eau", ylab="valeurs") #valeurs extrêmes avec des valeurs < 100 environs 5
boxplot(APICleaned$Fibres.alimentaires..g.100.g., main="Fibre alimentaires")
hist(APICleaned$Eau..g.100.g., )


# ==================================================================== #
# ************************ PART 2 *********************************** #
# ********************** ANALYSE BIDIMENSIONNELLE ****************** #
# =================================================================== #

pairs(APICleaned, c(APICleaned$Eau..g.100.g., APICleaned$Protéines..N.x.facteur.de.Jones..g.100.g., APICleaned$Protéines..N.x.6.25..g.100.g.,
                    APICleaned$Glucides..g.100.g., APICleaned$Lipides..g.100.g., APICleaned$Fibres.alimentaires..g.100.g., APICleaned$Polyols.totaux..g.100.g.,
                    APICleaned$Cendres..g.100.g., APICleaned$Alcool..g.100.g., APICleaned$Acides.organiques..g.100.g., APICleaned$AG.saturés..g.100.g.,
                    APICleaned$AG.monoinsaturés..g.100.g., APICleaned$AG.polyinsaturés..g.100.g., APICleaned$Sel.chlorure.de.sodium..g.100.g., APICleaned$Sodium..mg.100.g.))

#quoi de mieux pour une analyse bidimensionnelle qu'un scatterplot via la fonction pairs de R
pairs(APICleaned, pch = 21, bg = c("red", "green", "blue"))


# ==================================================================== #
# ************************ PART 3 *********************************** #
# ********************** DESCRIPTION MULTIVARIEE: ACP ****************** #
# =================================================================== #

#Center and Reduct
CenterReduct = function(Mat){
  return (scale(Mat, center = TRUE, scale=TRUE)*sqrt(nrow(Mat)/(nrow(Mat)-1)))
}

APICleaned =  subset(APICleaned, select = -c(Alcool..g.100.g.)) #delete of alcool column it provoques a NAN values
cleanData = subset(cleanData, select = -c(Alcool..g.100.g.)) #Data with the mean chease

dim(cleanData)

dataMatrix = data.matrix(APICleaned) #transform APICleaned to Matrix
dataMatrix2 = data.matrix(cleanData)

Inertia = function(Mat){
  Q=matrix(0,nrow=ncol(Mat),ncol=ncol(Mat))
  diag(Q)=1
  D=matrix(0,nrow=nrow(Mat),ncol=nrow(Mat))
  diag(D)=1/nrow(Mat)
  return(t(Mat)%*%D%*%Mat%*%Q)
}


#valeurs propres
EigenValues = function(Mat){ 
  return(eigen(Mat)$values)
}

#vecteurs propres
EigenVectors = function(Mat){
  return(eigen(Mat)$vectors)
  
}

#Calcul des coordonnées des individus
Fcoordonate = function(Mat,vec){
  coord = matrix(0,nrow=nrow(Mat),ncol=ncol(vec))
  for(i in 1:ncol(vec)){
    coord[,i]=Mat%*%vec[,i]
  }
  return(coord)
}

#Calcul des coordonnées des variables
Gcoordonate = function(val,vec){
  coord=matrix(0,nrow=nrow(vec),ncol=length(val))
  for(i in 1:ncol(vec)){
    coord[,i]=sqrt(val[i])*vec[,i]
  }
  return(coord)
}

#contribution individu
indContribution = function(Mat, D, EigValues, Fco, axesNb){
  
  result = matrix(0, nrow = nrow(Mat), ncol = axesNb)
  for(i in 1:dim(Mat)[1]){
    for(k in 1:axesNb){
      result[i,k] = (D[i,i]*Fco[i,k]**2)/EigValues[k]
      #result[i,k] = (Fco[i,k]**2)/dim(Mat)[1]*EigValues[k]
    }
  }
  return(result)
}

#contribution variable
varContribution = function(Mat, Gco, EigValues, axesNb){
  result = matrix(0, nrow = dim(Mat)[2], ncol = axesNb)
  for(j in 1:dim(Mat)[2]){
    for(k in 1:axesNb){
      result[j,k] = (Gco[j,k]**2)/EigValues[k]
      #result[j,k] = (Gco[j,k]**2)/EigValues[k]
    }
  }
  return(result)
}


#calcul qualité representation des individus
indRepresentationQuality = function(Mat, Fco, axesNb){
  result = matrix(0, nrow = dim(Mat)[1], ncol = axesNb)
  for(i in 1:dim(Mat)[1]){
    somme = sum(Fco[i,]**2)
    
    for(k in 1:axesNb){
      result[i,k] = (Fco[i,k]**2)/somme
    }
  }
  return(result)
}

#calcul qualité representation des variables
varRepresentationQuality = function(Mat, Gco, axesNb){
  result = matrix(0, nrow = dim(Mat)[2], ncol = axesNb)
  for(j in 1:dim(Mat)[2]){
    somme = sum(Gco[j,]**2)
    
    for(k in 1:axesNb){
      result[j,k] = (Gco[j,k]**2)/somme
    }
  }
  return(result)
}


library("ade4")
ACP = function(Mat){
  
  Q=matrix(0,nrow=ncol(Mat),ncol=ncol(Mat)) #Qp Matrix
  diag(Q) = 1
  
  D=matrix(0,nrow=nrow(Mat),ncol=nrow(Mat)) #D = 1/n * In Matrix
  diag(D) = 1/nrow(Mat)
  
  
  Mcr = CenterReduct(Mat)
  View(Mcr)
  Mi = Inertia(Mcr)
  View(Mi)
  Ig = sum(diag(Mi)) #trace de la matrice d'inertie = trace matrice de corrélation = p
  
  EigValues = EigenValues(Mi)
  View(EigValues)
  EigValues
  barplot(EigValues)  
  plot(EigValues)
  
  EigVects = EigenVectors(Mi)
  View(EigVects)
  
  print("Pourcentages cumulées en fonction du nombre d'axe")
  print(round(cumsum((EigValues/sum(diag(Mi))) * 100), 2))
  
  axes = readline("Combien d'axes souhaitez-vous gardé ? ")
  axes = as.numeric(axes)
  
  if( !(1 <= axes && axes <= dim(Mcr)[2]) ){
    print("Le nombre d'axes ne peut être inférieur à 1 ou supérieur au nombre de colonnes")
    #print("Le nombre d'axe est fixé de ce pas à 2 par défaut")
    
  }else{
    #Coordonnées Individu
    Fco = Fcoordonate(Mcr, EigVects)
    View(Fco)
    
    #Coordonnées Variables
    Gco = Gcoordonate(EigValues, EigVects)
    View(Gco)
    
    #Contribution individus
    indContrib = indContribution(Mat, D, EigValues, Fco, axes)
    View(indContrib)
    
    #Contribution variables
    varContrib  = varContribution(Mat, Gco, EigVects, axes)
    View(varContrib)
    
    #Qualité de la representation for les individus
    indQuality = indRepresentationQuality(Mat, Fco, axes)
    View(indQuality)
    
    #Qualité de la representaiton pour les variables
    varQuality = varRepresentationQuality(Mat, Gco, axes)
    View(varQuality)
    
    
    ### Calcul de la somme des contributions des individus pour chaque axe retenu
    
    for(i in 1:axes){
      message("La somme des contributions des individus sur l'axe",i)
      print(sum(indContrib[,i]))
    }
    
    #Graphique
    #Visualisation des individus
    s.label(Fco, xax=1, yax=2, label = 1:dim(Fco)[1], clabel = 0.7)
    
    #Visualisation des variables
    s.corcircle(Gco, xax =1, yax=2, label = 1:dim(Gco)[1], clabel = 0.7)
    
    # plot(Fco,xlab="F1",ylab="F2")
    # text(Fco,cex=0.65,pos=3,labels=1:nrow(Mat))
    
    
  }
  
  
}

ACP(dataMatrix)

ACP(dataMatrix2) #ACP after adding the average food to the dataset


