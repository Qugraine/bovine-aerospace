#Test scripts for furiousblue package

#Make sure balloon physics are sensible
source("../furiousblue/R/FlighPhysics.R")

#Test a hydrogen balloon at sea level

#Balloon specs
mb <- 1 #balloon weight (kg)
mi <- 1 #package weight (kg)
L <- 1.5 #Lift measured in the field
po <- 1 #Overpressure coefficient
ma <- 0.02897 #Molar mass of air
mp <- 1.008 * 2/1000 #Molar mass of lift gas (hydrogen, in this case)
max.radius <- 3 #Maximum allowed radius (m)
sealed <- TRUE #Not open to atmosphere
T.diff <- 0 #Temperature difference between balloon and atmosphere
acd <- 0.5 #Coefficient of drag
T.interior <- NULL #Internal temperature

#Atmosphere specs
pa <- 100000 #Pressure (Pascals)
T.ambient <- 293 #Temperature (K)

Mp <- MolsFromLift(mb, L, po, ma, mp) #Get mols of gas in balloon
v <- AscentVelocity(mi, Mp, mp, ma, pa, po, acd, max.radius, T.ambient, T.diff = T.diff, T.interior = T.interior, sealed = sealed)
