#Code for calculating balloon flight and payload descent

AscentCd <- function(Re = NULL, formula = "gallice") {
    #Determine the coefficient of drag of an ascending balloon.
    #This can be done in many ways, some of which take Reynold's number into account.
    #INPUTS
    #    RE is Reynold's number
    #    formula is the formula to use to calculate the drag coefficient (see notes below)
    #OUTPUTS
    #    CD is the coefficient of drag    

    #smoothsphere taken from http://www.engineeringtoolbox.com/drag-coefficient-d_627.html

    #gallice taken from Gallice et al (2011) "Modeling the ascent of sounding balloons: derivation of the vertical air motion" Atmospheric Measuring Techniques (4) 2235-2253
  #They note it is only strictly valid for a certain type of weather balloon at night

    #mikhailov taken from Mikhailov and Freire (2013) "The drag coefficient of a sphere: An approximation using Shanks transform", Powder Technology 237 432-435

    Cd <- NULL
    if(formula == "smoothsphere") { #From Engineering Toolbox website
        Cd <- 0.5
    }

    if(formula == "gallice")  {
        Cd <- 0.04808 * (log(Re))^(2) - 1.406 * log(Re) + 10.49
    }

    if(formula == "mikhailov") {
        Cd <- (777/(646 * Re)) * ((669806/875) + (114976/1155) * Re + (707/1380)*Re^2)/((32869/952) + (924/643) * Re + (1/385718) * Re^2)
    }

    if(formula == "smoothsphere" | is.null(Cd)) { #From Engineering Toolbox website
        Cd <- 0.5
    }

    invisible(Cd)
}

AscentRe <- function(ma, pa, Vb, T, va, r) {
    #Calculate Reynold's number for ascending balloon
    #INPUTS
    #   MA is molar mass of air (kg) 
    #   PA is air pressure (pa)
    #   T is temperature (K)
    #   VB is the volume of the balloon (m3)
    #   VA is ascent velocity (m/s)
    #   R is balloon radius (m)
    #OUTPUTS
    #   RE is the Reynold's number of the balloon

    #Calculate dynamic viscosity of air using Sutherland's formula

    mu0 <- 0.00001827 #Reference viscosity (pa s)
    To <- 291.15 #Reference temperature (K)
    C <- 120 #Sutherland's constant for air
    R <- 8.3145 #Universal gas constant

    dm <- (pa/(R*T))*ma #Air density

    mu <- mu0 * ((To + C)/(T + C)) * (T/To)^(1.5)

    invisible(dm*va*r*2.0/mu) #Reynold's number
}

AscentVelocity <- function(mi, Mp, mp, ma, pa, po, acd, max.radius, T.ambient, T.diff = 0, T.interior = NULL, sealed = TRUE) {
    #Calculate the ascent velocity of the balloon.
    #INPUTS
    #    MI is mass of instrument package, parachute, etc (kg)
    #    Mp is mols of lift gas
    #    MP is molar mass of lift gas (kg)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    PA is ambient air pressure (pascals)
    #    PO is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
    #    ACD is coefficient of drag (unitless)
    #    MAX.RADIUS is the radius of a non elastic envelope or the burst radius of an elastic envelope
    #    T.AMBIENT is temperature outside of balloon (kelvin)
    #    T.DIFF is temperature difference between balloon and outside (kelvin)
    #    T.INTERIOR is absolute interior temperature (overrides T.DIFF) (kelvin)
    #    SEALED denotes if the envelope is open to the atmosphere (and thus can lose gas) or is sealed shut (gas is retained)
    #OUTPUTS
    #    VA is ascent velocity (m/s)

    if(T.diff != 0 & !is.null(T.interior)) { #Note - this should be returned only once, not on every iteration
         warning("A temperature difference from ambient and a fixed internal temperature are both defined for the balloon.  The fixed interior temperature shall override the temperature difference.")
    }

    R <- 8.3145 #Universal gas constant
    g <- 9.81 #Gravitational acceleration

    if(is.null(T.interior)) {
        rv <- BalloonSphere(Mp, pa, po, T.ambient + T.diff)
    } else {
        rv <- BalloonSphere(Mp, pa, po, T.interior)
    }

    if(max.radius < rv$r) {
        rv$r <- max.radius
        rv$V <- (4/3)* pi * max.radius ^ 3
    }

    disp.mass=((pa*rv$V)/(R*T.ambient))*ma #Mass of displaced fluid

    if(sealed) {
        prop.mass=Mp*mp #Mass of lift gas
    } else {
        if(is.null(T.interior)) {
            prop.mass <- ((pa*rv$V)/(R*(T.ambient + T.diff)))*mp
        } else {
            prop.mass <- ((pa*rv$V)/(R*(T.interior)))*mp
        }
    }

    L <- (disp.mass-prop.mass-mi)*g #Net lift force

    A <- pi*rv$r^2 #Cross sectional area of balloon

    dm <- disp.mass/rv$V #Air density
    va <- sqrt((2*L)/(acd*dm*A))
    invisible(va) #Ascent velocity
}

BalloonSphere <- function(Mp, pa, po, T) {#Get balloon dimensions
   #Calculate radius and volume of balloon if it is elastic
   #INPUTS
   #   MP is mols of lift gas
   #   PA is ambient air pressure (pascals)
   #   po is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
   #   T is internal temperature (kelvin)
   #OUTPUTS
   #   A list with elementds
   #   r - Balloon radius (m)
   #   V - Balloon volume (m3)

   R <- 8.3145 #Universal gas constant

   V <- (Mp*R*T)/(pa*po) #Balloon volume (m3)

    numerator <- 3.0*Mp*R*T

    denominator <- 4*pi*(pa*po)



    r <- ((3 * V)/(4 * pi)) ^(1/3)
    invisible(list(r = r, V = V))
}

DescentVelocity <- function(mi, ma, pa, rp, T, dcd) {
    #Calculate velocity of a (non spherical) descending object
    #INPUTS
    #    MI is mass of instrument package, parachute, etc (kg)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    PA is ambient air pressure (pascals)
    #    RP is parachute radius (assume circular) (m)
    #    T is temperature (Kelvin)
    #    DCD is coefficient of drag (unitless)

    #CODE VERIFIED BY NASA ROCKET PARACHUTE DESCENT MODEL

    R <- 8.3145 #Universal gas constant
    g <- 9.81 #Gravitational acceleration

    dm <- pa*ma/(R*T) #Air density
    A <- pi*rp*rp #Cross sectional area of parachute

    vd <- sqrt(2*mi*g/(dcd*dm*A)) #Descent velocity
    invisible(vd)
}

MolsFromLift <- function(mb, L, po, ma, mp) {
    #Calculate how many mols of lift gas you added to your balloon to get the measured lift at the ground
    #This assumes ISOTHERMAL CONDITIONS...so the balloon is the same temperature as the ambient air
    #INPUTS
    #    MB is the mass of the balloon envelope (kg)
    #    L is the amount of lift measured in the field (kg)
    #    PO is the overpressure (how much higher the pressure is inside the balloon, as a dimensionless ratio balloon pressure/ambient air pressure)
    #    MA is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #    MP is the molar mass of the propellant
    #OUTPUTS
    #MP  how many mols of lift gas you have added to get measured lift

    total.mass <- L + mb

    Mp <- total.mass/((ma/po)-mp)

    invisible(Mp) #How many mols of lift gas you added to your balloon to get the measured lift
}

MolsFromVolume <- function(V, pa, T, T.diff = 0, T.interior = NULL) {
   #How many mols of gas are inside the balloon, assuming the balloon is open to the atmosphere
   #For example, how many mols of air are inside my fully inflated solar balloon

   if(T.diff != 0 & !is.null(T.interior)) {
         warning("A temperature difference from ambient and a fixed internal temperature are both defined for the balloon.  The fixed interior temperature shall override the temperature difference.")
    }

   R <- 8.3145 #Universal gas constant

   if(is.null(T.interior)) {
       Mp <- pa * V / (R * (T + T.diff))
   } else {
       Mp <- pa * V / (R * T.interior)
   }

   invisible(Mp)
}

