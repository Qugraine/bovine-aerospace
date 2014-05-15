launch_lat=35.823833
launch_lon=-79.105557
model_parameters=[10000, 0.1, 0.1, 0.1, 0.1, 30000, 5000, 40000]
#model_parameters[0]-Number of times to run the monte carlo model
#model_parameters[1]-standard deviation of total ascent velocity
#model_parameters[2]-standard deviation of total descent velocity
#model_parameters[3]-standard deviation of x windspeed
#model_parameters[4]-standard deviation of y windspeed
#model_parameters[5]-Average burst height
#model_parameters[6]-Burst height standard deviation
#model_parameters[7]-Maximum ascent elevation allowed

def ascent_reynolds_number(ma, pa, Vb, T, va, r):
    #Calculate Reynold's number of an ascending sphere
    #ma is mass of air (kg)
    #pa is air pressure (pa)
    #Vb is the volume of the balloon (m3)
    #T is temperature (K)
    #va is ascent velocity (m/s)
    #r is balloon radius

    import math

    R=8.3145 #Universal gas constant
    uo=0.00001827 #Reference viscosity (pa s)
    To=291.15 #Reference temperature (K)
    C=120 #Sutherland's constant for air

    dm=(pa/(R*T))*ma #Air density
    #u=uo*((To+C)/float(T+C))*(math.pow(T/float(To),3/2.0)) #Dynamic viscosity (Sutherland's Formula)
    u=(uo*(To+C)*math.pow(T/To,1.5))/(T+C)

    Re=(dm*va*r*2.0)/u #Reynold's number

    return Re

def ascent_cd(Re):
    #Calculate coefficient of drag of an ascending sphere
    #Re=reynold's number

    import math
    acd=(24/Re)+((2.6*(Re/5.0))/(1+math.pow(Re/5.0, 1.52)))+(0.411*math.pow(Re/263000.0,-7.94)/(1+math.pow(Re/263000.0,-8.00)))+(math.pow(Re, 0.80)/461000.0)
    return acd

    

def calculate_lat_lon(start_lat, start_lon, x, y):
    #Start lat is the latitude of the point you started from
    #Start lon is the longitude of the point you started from
    #x is how many meters you have gone east-west (m)
    #y is how many meters you have gone north south (m)
    
    dlen=111325 #Length of a degree of latitude (m)
    import math
    lat2=start_lat+(y/dlen) 
    lon2=start_lon+(x/(dlen*(math.cos((math.pi/180)*(lat2+start_lat)/2))))
    return [lat2, lon2]

def get_most_recent_model_time():
    #This function grabs the most recent available model date from the University of Wyoming Balloon Model page.
    import urllib
    f=urllib.urlopen("http://weather.uwyo.edu/polar/balloon_traj.html")
    s=f.read()
    f.close()
    chopstr=s.split("option")
    chopstr=chopstr[1]
    chopstr=chopstr.split("\"")
    date=chopstr[3]
    return date

def mols_from_lift(mb, L, po, ma, mp):
    #mb is the mass of the balloon envelope (kg)
    #L is the amount of lift measured in the field (kg)
    #po is the overpressure (how much higher the pressure is inside the balloon, as a dimensionless ratio balloon pressure/ambient air pressure)
    #ma is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #mp is the molar mass of the propellant

    total_mass=L+mb

    Mp=float(total_mass)/((ma/float(po))-mp)

    return Mp #How many mols of lift gas you added to your balloon to get the measured lift

def balloon_sphere(Mp, pa, po, T):
    #Mp is mols of lift gas
    #pa is ambient air pressure (pascals)
    #po is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
    #T is temperature (kelvin)

    import math

    R=8.3145 #Universal gas constant
    V=float(Mp*R*T)/(pa*po) #Balloon volume (m3)
    numerator=3.0*Mp*R*T
    denominator=4*math.pi*(pa*po)

    r=math.pow(numerator/denominator, (float(1)/3)) #Balloon radius (m)
    return r, V
    

def ascent_velocity(mi, Mp, mp, ma, pa, po, T, acd):
    #mi is mass of instrument package, parachute, etc (kg)
    #Mp is mols of lift gas
    #mp is molar mass of lift gas (kg)
    #ma is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #pa is ambient air pressure (pascals)
    #po is balloon overpressure, ratio of balloon internal pressure to ambient air pressure
    #T is temperature (kelvin)
    #acd is coefficient of drag (unitless)

    #CONFIRMED VIA EMPIRICAL DATA

    import math

    R=8.3145 #Universal gas constant
    g=9.81 #Gravitational acceleration
    
    rb, Vb=balloon_sphere(Mp, pa, po, T) #Get balloon dimensions

    disp_mass=((pa*Vb)/(R*T))*ma #Mass of displaced fluid
    prop_mass=Mp*mp #Mass of lift gas
    L=(disp_mass-prop_mass-mi)*g #Net lift force
    A=math.pi*rb*rb #Cross sectional area of balloon
    dm=disp_mass/Vb #Air density

    va=math.pow((2*L)/(acd*dm*A), 0.5)
    return va #Ascent velocity
    
def descent_velocity(mi, ma, pa, rp, T, dcd):
    #mi is mass of instrument package, parachute, etc (kg)
    #ma is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #pa is ambient air pressure (pascals)
    #rp is parachute radius (assume circular) (m)
    #T is temperature (Kelvin)
    #dcd is coefficient of drag (unitless)

    #CODE VERIFIED BY NASA ROCKET PARACHUTE DESCENT MODEL
    
    import math

    R=8.3145 #Universal gas constant
    g=9.81 #Gravitational acceleration

    dm=float(pa)*ma/(R*T) #Air density
    A=math.pi*rp*rp #Cross sectional area of parachute

    vd=math.pow(2*mi*g/(dcd*dm*A), 0.5) #Descent velocity


    return vd
    
def get_trajectory_from_internet(start_lat, start_lon, max_altitude, model_start_time, forecast, model_format):
    #start_lat is the latitude from which the balloon was launched
    #start_lon is the longitude from which the balloon was launched
    #max_altitude is the altitude where the University of Wyoming assumes the balloon will burst (in meters)
    #model_start_time is the date/timestring of the GFS model yyyymmddhh (hours are 00, 06, 12, 18 in GMT)
    #forecast is how many hours ahead to project the model before data is extracted
    #Acceptable forecast inputs are 0 (latest model run), 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84 hours ahead
    #Input MODEL FORMAT is "list" to save the model as a text file and return an input for monte_carlo_balloon_flight
    #Input MODEL FORMAT as "kml" to download a KML file from the University of Wyoming website
    #
    #
    #Returns a 2 element vector: [lat, lon] of the launch point
    #Also returns a tuple containing the input for the monte_carlo_balloon_flight model

    #This function returns an elevation profile for windspeed.
    import pickle
    import urllib
    import time
    from datetime import datetime

    accepted_forecasts=[0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84]
    if forecast not in accepted_forecasts:
        print("Could not find a forecast for the requested time.  Acceptable forecasts are:")
        for hr in accepted_forecasts:
            print(str(hr)+" "),
        print("\n hours ahead of the most recent model run time")
        return 1, 1

    now=datetime.utcnow() #Get forecast for right now in GMT
    year=str(now.year)
    month=str(now.month)
    day=str(now.day)
    hour=str(now.hour)
    minute=str(now.minute)
    if len(month)==1:
        month="0"+month
    if len(day)==1:
        day="0"+day
    if len(hour)==1:
        hour="0"+hour
    if len(minute)==1:
        minute="0"+minute

    attempts=0 #How many times we tried to connect to model
    max_attempts=10 #How many times to try before giving up
    pause=30 #How many seconds to pause in between attempts
    urlstring="http://weather.uwyo.edu/cgi-bin/balloon_traj?TIME="+model_start_time+"&FCST="+str(forecast)+"&POINT=none&LAT="+str(start_lat)+"&LON="+str(start_lon)+"&TOP="+str(max_altitude)+"&OUTPUT="+model_format+"&Submit=Submit&.cgifields=POINT&.cgifields=FCST&.cgifields=TIME&.cgifields=OUTPUT"
    while 1:
        if attempts<max_attempts:
            attempts+=1
            print("Waiting for response from server...")
            f=urllib.urlopen(urlstring)
            print("Response received.  Analyzing response...")
            s=f.read()
            if "sorry" in s.lower():
                print("Server is too busy to respond; waiting "+str(pause)+" seconds and trying again...")
                time.sleep(pause) #Wait ten seconds, then try again
            else:
                if model_format=="list":
                    model_file=open(str(start_lat).replace(".","-")+str(start_lon).replace(".","-")+"_"+year+month+day+hour+".txt", "w")
                    model_file.write(s)
                    model_file.close()
                    print("Model file "+str(start_lat).replace(".","-")+str(start_lon).replace(".","-")+"_"+year+month+day+hour+".txt written as HTML-formatted text file!")
                    flight_data=read_uw_model(s)
                    break
                else:
                    model_file=open(str(start_lat).replace(".","-")+str(start_lon).replace(".","-")+"_"+year+month+day+hour+".kml", "w")
                    model_file.write(s)
                    model_file.close()
                    print("Model file " + str(start_lat).replace(".","-")+str(start_lon).replace(".","-")+"_"+year+month+day+hour+".kml written as KML file!")
                    flight_data=(0, 0, 0, 0, 0, 0, 0)
                    break
        else:
            print("Maximum attempt limit reached.  Try again some other time!")
            flight_data=(1, 1, 1, 1, 1, 1, 1)
            break
    return([start_lat, start_lon], flight_data)

def read_uw_model(f):
    #Read HTTP code from UW model to build a vertical forecast
    #f is the file containing weather data (either loaded from disk, or returned from the University of Wyoming website

    
    try:
        s=f.read()
    except:
        s=f
    sline=s.split("<pre>") #Get trajectory data from HTTP string
    sline=sline[2]
    sline=sline.split("</PRE>")
    sline=sline[0]
    trajectory=sline.split("\n")


    elev=[] #Elevation (m)
    x=[] #East-West Windspeed (m/s)
    y=[] #North-South windspeed (m/s)
    z=[] #Ascent/Descent velocity (m/s), will be depreciated
    t=[] #Elapsed time (s)
    p=[] #Pressure (pascals)
    T=[] #Temperature (kelvin)
    for step in trajectory:
        if step is not "":
            straj=step.split()
            flight_time=straj[0].split(':')
            t.append(int(flight_time[0])*60*60+int(flight_time[1])*60+int(flight_time[2]))
            x.append(float(straj[6]))
            y.append(float(straj[7]))
            z.append(float(straj[8]))
            p.append(float(straj[9])*100)
            T.append(float(straj[10])+273.15)
            elev.append(float(straj[3]))
    weather=(x, y, z, t, p, T, elev)    
    print("University of Wyoming model data successfully read!")
    try:
        f.close()
    except:
        pass
    return weather

def build_model(weather, monte_carlo):
    #Interpolate weather data from UW into a 2d model of Temperature (K), Pressure (pascals), East Windspeed (m/s), North Windspeed (m/s)
    #monte_carlo is a 4 element vector with standard deviation data to run a monte-carlo model [x, y, p, T] where x and y are E-W and N-S windspeed
    #adjusted by a fraction of the total magnitude of each (m/s); p is standard deviation in pressure (pascals), T is standard deviation in degrees K.
    #Model is in 100 m elevation increments from 0 to 50000

    import numpy as np
    from scipy import interpolate
    import random

    #import matplotlib.pyplot as plt

    #Get variables
    x=weather[0] #Windspeed x (m/s)
    y=weather[1] #Windspeed y (m/s)
    T=weather[5] #Temperature (K)
    p=weather[4] #Pressure (pa)
    elev=weather[6] #Elevation (m)
    max_elev_ind=elev.index(max(elev))+1
    #Only use the ascent portion of the model

    elev=elev[:max_elev_ind]
    x=x[:max_elev_ind]
    y=y[:max_elev_ind]
    T=T[:max_elev_ind]
    p=p[:max_elev_ind]

    for k in range(len(x)):
        x[k]=x[k]+x[k]*random.gauss(0, monte_carlo[0]) #calculate adjustment based on standard deviation factor as fraction of magnitude of x
        y[k]=y[k]+y[k]*random.gauss(0, monte_carlo[1]) #calculate adjustment based on standard deviation factor as fraction of magnitude of y
        p[k]=p[k]+random.gauss(0, monte_carlo[2]) #calculate adjustment based on standard deviation in pascals
        T[k]=T[k]+random.gauss(0, monte_carlo[3]) #Calculate adjustment based on standard deviation in degrees Kelvin

    xtck=interpolate.splrep(elev, x, s=0)
    #xint=interpolate.splev(np.arange(elev[0],elev[-1],100),xtck,der=0)
    ytck=interpolate.splrep(elev, y, s=0)
    #yint=interpolate.splev(np.arange(elev[0],elev[-1],100),ytck,der=0)
    Ttck=interpolate.splrep(elev[:max_elev_ind], T, s=0)
    #Tint=interpolate.splev(np.arange(elev[0],elev[-1],100),Ttck,der=0)
    ptck=interpolate.splrep(elev[:max_elev_ind], p, s=0)
    #pint=interpolate.splev(np.arange(elev[0],elev[-1],100),ptck,der=0)    
    #plt.plot(xint, np.arange(elev[0],elev[-1],100))
    #plt.plot(x, elev, 'ro')
    #plt.show()

    return (xtck, ytck, Ttck, ptck, elev)

def launch(launch_point, launch_time, ascent_model, descent_model, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, file_name):
    #Runs balloon trajectory model
    #launch point is where the balloon is launched from ([lat lon])
    #launch_time is when the balloon is launched (currently GMT)
    #ascent model is the forecast-based interpolation of windspeed (m/s), temperature (K), and pressure (pa) from up to 30000-40000 m at launch point
    #ascent model is the forecast-based interpolation of windspeed (m/s), temperature (K), and pressure (pa) from up to 30000-40000 m at expected landing point
    #model_resolution is the resolution of the vertical model (m)
    #forward determines if the model should be run forward (landing point given launch point) or reverse (launch point given landing point)
    #mb is the balloon mass (kg)
    #L is the field-measured lift prior to attaching instrument package (kg)
    #po is the ratio of balloon overpressure to ambient air (unitless), describes how much higher the pressure is inside of the balloon vs the outside
    #ma is molar mass of displaced fluid (air, in this case, if dry it is 0.02897) (kg)
    #mp is the molar mass of lift gas (kg)
    #mi is the mass of the instrument package, parachute, etc (not including balloon mass) (kg)
    #dcd is the coefficient of drag of a descending instrument package, hopefully with parachute! (dimensionless)
    #rp is the radius of the parachute (m)
    #rb is the burst radius of the balloon (how big the balloon gets just before bursting) (m)
    #file_name is the name of the data text file generated by this function

    #Import necessary libraries

    from scipy import interpolate as i
    #from matplotlib import pyplot as plt
    import math

    #Open output file, write model parameters
    
    out=open(file_name, "w")
    out.write("BALLOON FLIGHT MODEL VERSION 1.0\n")
    out.write("BOVINE AEROSPACE in association with CHICKEN HURLERS INCORPORATED\n")
    out.write("Start Point (lat lon dec degree):\t"+str(launch_point[0])+"\t"+str(launch_point[1])+"\n")
    out.write("Start Time: (UTM)\t"+str(launch_time)+"\n")
    out.write("Model vertical resolution (m):\t"+str(model_resolution)+"\n")
    if forward:
        out.write("Model direction:\tFORWARD (launch to landing)\n")
    else:
        out.write("Model direction:\tREVERSE (landing to launch)\n")
    out.write("Balloon mass (kg):\t"+str(mb)+"\n")
    out.write("Measured lift (before instrument package is attached) (kg):\t"+str(L)+"\n")
    out.write("Balloon overpressure (ratio of interior to exterior pressure) (unitless):\t"+str(po)+"\n")
    out.write("Molar mass of air (kg):\t"+str(ma)+"\n")
    out.write("Molar mass of lift gas (kg):\t"+str(mp)+"\n")
    out.write("Mass of instrument package and parachute (kg):\t"+str(mi)+"\n")
    out.write("Coefficient of drag of descending instrument package and parachute (unitless):\t"+str(dcd)+"\n")
    out.write("Parachute radius (m):\t"+str(rp)+"\n")
    out.write("Balloon burst radius (m)\t"+str(rb)+"\n")

    #Determine molar mass given initial lift

    Mp=mols_from_lift(mb, L, po, ma, mp)

    #Ascend

    xtck=ascent_model[0]
    ytck=ascent_model[1]
    Ttck=ascent_model[2]
    ptck=ascent_model[3]
    elev=ascent_model[4]
    
    x_pos=0 #Distance traveled (East-West) (m)
    y_pos=0 #Distance traveled (North-South) (m)
    z_pos=elev[0] #Elevation (m)
    t_pos=0 #Elapsed time (s)
    latlon_position=launch_point

    xy=[] #Final distance traveled, in X (East-West) and Y (North South)


    cs=model_resolution/2
    layers=0
    av_bin=[]
    alev_bin=[]
    out.write("Lat\tLon\tElev (m)\tTime (s)\t East-West Windspeed (m\\s)\tNorth-South Windspeed (m\\s)\tAscent Velocity (m\\s)\tTemp (k)\t Pressure (pa)\tBalloon Radius (m)\n")
    out.write("MODEL START\n")
    out.write("ASCEND\n")
    Re_bin=[]
    cd_bin=[]
    x_bin=[]
    y_bin=[]
    
    T=i.splev(elev[0]+1,Ttck,der=0)
    pa=i.splev(elev[0]+1, ptck, der=0)
    acd=0.4 #Starting ascent coefficient of drag
    #Iterate to starting ascent coefficient of drag
    for k in range(25):
        r, V=balloon_sphere(Mp, pa, 1.1, T)
        av=ascent_velocity(mi, Mp, mp, ma, pa, po, T, acd)
        Re=ascent_reynolds_number(ma, pa, V, T, av, r)
        acd=ascent_cd(Re)

    #Calculate ascent trajectory
    while 1:
        z_pos+=cs #Calculate midpoint of layer for p, T, x, and y
        T=i.splev(z_pos,Ttck,der=0)
        pa=i.splev(z_pos, ptck, der=0)
        x=i.splev(z_pos, xtck, der=0)
        x_bin.append(x)
        y=i.splev(z_pos, ytck, der=0)
        y_bin.append(y)
        r, V=balloon_sphere(Mp, pa, 1.1, T) #Calculate the balloon radius at the midpoint of the layer
        av=ascent_velocity(mi, Mp, mp, ma, pa, po, T, acd) #Calculate the ascent rate at the midpoint of the layer
        av_bin.append(av)
        Re=ascent_reynolds_number(ma, pa, V, T, av, r) #Calculate reynold's number
        acd=ascent_cd(Re)
        Re_bin.append(Re)
        cd_bin.append(acd)
        t=model_resolution/float(av) #Calculate the amount of time spent in each layer
        t_pos+=t
        x_pos+=x*t #How far in the x direction the balloon travels in the layer
        y_pos+=y*t #How far in the y direction the balloon travels in the layer
        alev_bin.append(z_pos)        
        if r>rb: #if the radius of the balloon exceeds the burst radius, stop ascending
            break
        #print("Elevation: "+str(z_pos)+" p: "+str(math.sqrt(x_pos*x_pos+y_pos*y_pos)))
        z_pos+=cs
        layers+=1
        latlon_position=calculate_lat_lon(latlon_position[0], latlon_position[1], x*t, y*t)
        out.write(str(latlon_position[0])+"\t"+str(latlon_position[1])+"\t"+str(z_pos)+"\t"+str(t_pos)+"\t"+str(x)+"\t"+str(y)+"\t"+str(av)+"\t"+str(T)+"\t"+str(pa)+"\t"+str(r)+"\n")        
        if layers>295:
            print("layerbreak at elevation "+str(z_pos))
            break
    ascent_time=t_pos

    #Descend
        
    xtck=descent_model[0]
    ytck=descent_model[1]
    Ttck=descent_model[2]
    ptck=descent_model[3]
    elev=descent_model[4]

    magnitude=[]
    dv_bin=[]
    elev_bin=[]
    out.write("DESCEND\n")
    out.write("Lat\tLon\tElev (m)\tTime (s)\t East-West Windspeed (m\\s)\tNorth-South Windspeed (m\\s)\tDescent Velocity (m\\s)\tTemp (k)\t Pressure (pa)\tParachute Radius (m)\n")
    while 1:
        z_pos-=cs #Calculate midpoint of layer for p, T, x, and y
        T=i.splev(z_pos,Ttck,der=0)
        pa=i.splev(z_pos, ptck, der=0)
        x=i.splev(z_pos, xtck, der=0)
        y=i.splev(z_pos, ytck, der=0)
        #print(str(mi)+"\t"+str(rp)+"\t"+str(dcd))
        dv=descent_velocity(mi, ma, pa, rp, T, dcd) #Calculate the descent rate at the midpoint of the layer
        dv_bin.append(dv)
        elev_bin.append(z_pos)
        t=model_resolution/float(dv) #Calculate the amount of time spent in each layer
        t_pos+=t
        x_pos+=x*t #How far in the x direction the balloon travels in the layer
        y_pos+=y*t #How far in the y direction the balloon travels in the layer
        magnitude.append(math.sqrt(x*x+y*y))
        z_pos-=cs
        latlon_position=calculate_lat_lon(latlon_position[0], latlon_position[1], x*t, y*t)
        out.write(str(latlon_position[0])+"\t"+str(latlon_position[1])+"\t"+str(z_pos)+"\t"+str(t_pos)+"\t"+str(x)+"\t"+str(y)+"\t"+str(dv)+"\t"+str(T)+"\t"+str(pa)+"\t"+str(rp)+"\n")  
        if z_pos<=elev[1]:
            break
    out.write("MODEL END")
    landing_point=calculate_lat_lon(launch_point[0], launch_point[1], x_pos, y_pos)
    if 1==0:
        print("Landing lat "+str(landing_point[0]))
        print("Landing lon "+str(landing_point[1]))
        print("Landing velocity "+str(dv_bin[-1])+" m/s ("+str(dv_bin[-1]*3.28*3600/5280)+" mph)")
        print("Maximum windspeed (mph) "+str(max(magnitude)*3.28*3600/5280))
        print("Ascent time (dec. hour) "+str(ascent_time/3600))
        print("Descent time (dec. hour) "+str((t_pos-ascent_time)/3600))
        print("Mols of helium in balloon "+str(Mp))
        stp_vol=(Mp*8.3145*293.15)/float(100000)
        print("Volume of helium at STP "+str(stp_vol)+" m3 ("+str(stp_vol*35.3145)+ " ft3)")
    out.close()
    xy.append(x_pos)
    xy.append(y_pos)
    #plt.plot(x_bin, alev_bin)
    #plt.plot(y_bin, alev_bin)
    #plt.show()
    return landing_point, xy

def build_flight_kml(file_name):
    #Write entire jake model to KML
    kml_file=open(file_name[:-4]+".kml","w")
    data_file=open(file_name, "r")
    text=data_file.read()
    flight_text=text.split("ASCEND")
    flight_text=flight_text[1]
    flight_text=flight_text.split("DESCEND")
    ascent_text=flight_text[0]
    descent_text=flight_text[1]
    kml_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    kml_file.write("<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n")
    kml_file.write("<Document>")
    kml_file.write("	<StyleMap id=\"msn_open-diamond\">\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>normal</key>\n")
    kml_file.write("			<styleUrl>#sn_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>highlight</key>\n")
    kml_file.write("			<styleUrl>#sh_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("	</StyleMap>\n")
    kml_file.write("	<Style id=\"sh_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.4</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")
    kml_file.write("	<Style id=\"sn_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.2</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")

    #loop over lines in file
    ascent_line=ascent_text.split("\n")
    ascent_line=ascent_line[1:-2]
    for ascent in ascent_line:
        ascent_array=ascent.split("\t")
        kml_file.write("<Placemark>\n")
        kml_file.write("<Name>Jake</Name>\n")
        kml_file.write("<visibility>1</visibility>\n")
        kml_file.write("<styleUrl>#msn_open-diamond</styleUrl>\n")
        kml_file.write("<Point>\n")
        kml_file.write("<altitudeMode>absolute</altitudeMode>\n")
        kml_file.write("<coordinates>"+str(ascent_array[1])+","+str(ascent_array[0])+","+ascent_array[2]+"</coordinates>\n")
        kml_file.write("</Point>\n")
        kml_file.write("</Placemark>\n")

    descent_line=descent_text.split("\n")
    descent_line=descent_line[1:-2]
    for descent in descent_line:
        descent_array=descent.split("\t")
        kml_file.write("<Placemark>\n")
        kml_file.write("<Name>Jake</Name>\n")
        kml_file.write("<visibility>1</visibility>\n")
        kml_file.write("<styleUrl>#msn_open-diamond</styleUrl>\n")
        kml_file.write("<Point>\n")
        kml_file.write("<altitudeMode>absolute</altitudeMode>\n")
        kml_file.write("<coordinates>"+str(descent_array[1])+","+str(descent_array[0])+","+descent_array[2]+"</coordinates>\n")
        kml_file.write("</Point>\n")
        kml_file.write("</Placemark>\n")
    kml_file.write("</Document>\n")
    kml_file.write("</kml>")
    kml_file.close()
    data_file.close()
    return 0

def mc2(launch_point, launch_time, weather, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, trials):
    #Run the jake model as a Monte Carlo simulation and write to KML

    import random

    x=[]
    y=[]

    x_sd=0.1 #Standard deviation in X
    y_sd=0.1 #Standard deviation in Y
    p_sd=0 #Standard deviation in pressure
    T_sd=2 #Standard deviation in temperature

    kml_file=open("mc2.kml","w")
    kml_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    kml_file.write("<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n")
    kml_file.write("<Document>")
    kml_file.write("	<StyleMap id=\"msn_open-diamond\">\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>normal</key>\n")
    kml_file.write("			<styleUrl>#sn_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>highlight</key>\n")
    kml_file.write("			<styleUrl>#sh_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("	</StyleMap>\n")
    kml_file.write("	<Style id=\"sh_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.4</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")
    kml_file.write("	<Style id=\"sn_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.2</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")
    for k in range(trials):
        ascent_model=build_model(weather, [x_sd, y_sd, p_sd, T_sd])
        descent_model=ascent_model
        #L=L+random.uniform(-0.3, 0.3)
        po=random.uniform(1.0, 1.2)
        dcd=dcd+random.uniform(-0.5, 0.5)
        #mi=mi+random.uniform(-0.3, 0.3)
        rp=rp+random.uniform(-0.1, 0.1)
        rb=rb+random.uniform(-0.3, 0.3)
        file_name=r"C:\bovine_aerospace\monte_carlo_trial_trajectories\windspeed_"+str(k)+".txt"
        landing_point, xy=launch(launch_point, launch_time, ascent_model, descent_model, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, file_name)
        build_flight_kml(file_name)
        kml_file.write("<Placemark>\n")
        kml_file.write("<Name>Jake</Name>\n")
        kml_file.write("<visibility>1</visibility>\n")
        kml_file.write("<styleUrl>#msn_open-diamond</styleUrl>\n")
        kml_file.write("<Point>\n")
        kml_file.write("<altitudeMode>clampToGround</altitudeMode>\n")
        kml_file.write("<coordinates>"+str(landing_point[1])+","+str(landing_point[0])+",0</coordinates>\n")
        kml_file.write("</Point>\n")
        kml_file.write("</Placemark>\n")
        x.append(xy[0]/1000.0)
        y.append(xy[1]/1000.0)
        if k%100==0:
            print(". "),
    kml_file.write("</Document>\n")
    kml_file.write("</kml>")    
    kml_file.close()

def monte_carlo_balloon_flight(launch_point, flight_data, model_parameters):
    #This script traces the balloon flight solely based on data from UW model, including ascent and descent velocity.
    #Attempt to vary everything (winds, ascent, descent, burst height)
    import pickle
    import math
    import random

    x_pos=0
    y_pos=0
    elapsed_time=0

    tests=model_parameters[0] #How many times to run the simulation
    ascent_velocity_sd=model_parameters[1] #Standard deviation of ascent velocity as a fraction of model layer ascent velocity (so when ascent_velocity_sd=0.1, the standard deviation of ascent velocity at each time step is 10% of the ascent velocity
    descent_velocity_sd=model_parameters[2] ##Standard deviation of descent velocity as a fraction of model layer ascent velocity
    x_sd=model_parameters[3] #Standard deviation fraction of x-directional wind velocity
    y_sd=model_parameters[4] #Standard deviation fraction of y-directional wind velocity
    expected_max_elev=model_parameters[5] #Where we expect the balloon to burst
    max_elev_sd=model_parameters[6] #Standard deviation of burst elevation
    elev_ceiling=model_parameters[7] #Above which the balloon will not go.

    lat1=launch_point[0] #Launch latitude
    lon1= launch_point[1]#Launch longitude
    t=flight_data[3] #Time since launch
    x=flight_data[0] #X Velocity
    y=flight_data[1] #Y Velocity
    z=flight_data[2] #Z Velocity
    p=flight_data[4]
    T=flight_data[5]
    elev=flight_data[6] #Elevation


    final_x=[] #Final X Position
    final_y=[] #Final Y Position
    timestep=0 #How long the balloon spends in a particular wind regime

    #Initiate KML

    kml_file=open("monte_carlo.kml","w")
    kml_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    kml_file.write("<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n")
    kml_file.write("<Document>")
    kml_file.write("	<StyleMap id=\"msn_open-diamond\">\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>normal</key>\n")
    kml_file.write("			<styleUrl>#sn_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("		<Pair>\n")
    kml_file.write("			<key>highlight</key>\n")
    kml_file.write("			<styleUrl>#sh_open-diamond</styleUrl>\n")
    kml_file.write("		</Pair>\n")
    kml_file.write("	</StyleMap>\n")
    kml_file.write("	<Style id=\"sh_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.4</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")
    kml_file.write("	<Style id=\"sn_open-diamond\">\n")
    kml_file.write("		<IconStyle>\n")
    kml_file.write("			<scale>1.2</scale>\n")
    kml_file.write("			<Icon>\n")
    kml_file.write("				<href>http://maps.google.com/mapfiles/kml/shapes/open-diamond.png</href>\n")
    kml_file.write("			</Icon>\n")
    kml_file.write("		</IconStyle>\n")
    kml_file.write("		<ListStyle>\n")
    kml_file.write("		</ListStyle>\n")
    kml_file.write("	</Style>\n")

    highest_windspeed=[]
    max_elevs=[]
    for j in range(tests):
        x_pos=0
        y_pos=0
        elapsed_time=0
        prev_elev=0
        elev.append(elev[0]) #Where the balloon starts
        ascent_velocity=random.gauss(0, ascent_velocity_sd) #Ascent velocity variation
        descent_velocity=random.gauss(0, descent_velocity_sd) #Descent velocity variation
        x_s=random.gauss(0, x_sd) #Systematic deviation in x
        y_velocity=random.gauss(0, y_sd) #Systematic deviation in y
        max_elev=random.gauss(expected_max_elev, max_elev_sd) #Variation in maximum elevation
        highest_windspeed.append(0)
        f=0
        if max_elev>elev_ceiling:
            max_elev=elev_ceiling
        max_elevs.append(max_elev)
        for k in range(len(t)):
            if elev[k]>max_elev:
                while elev[f]>max_elev:
                    f=f+1
                    
            if z[k+f]<0:
                timestep=abs((elev[k+1+f]-elev[k+f])/(z[k+f]+z[k+f]*descent_velocity))
               # timestep=abs((elev[k+1+f]-elev[k+f])/descent_velocity(p[k+f-1],T[k+f-1],0.7,3,5))
                #print(z[k+f]/terminal_velocity(p[k+f-1],T[k+f-1],0.7,3,5))
            else:
                timestep=abs((elev[k+f]-elev[k+f-1])/(z[k+f-1]+z[k+f-1]*ascent_velocity))
            if 1==0:
                x_pos=x_pos+(float(x[k+f])+(float(x[k+f])*x_velocity))*(float(timestep)) #Systamatic wind variation throughout launch
                y_pos=y_pos+(float(y[k+f])+(float(y[k+f])*y_velocity))*(float(timestep)) #Systematic wind variation throughout launch
            else:
                x_velocity=random.gauss(float(x[k+f]), float(x[k+f])*x_sd)
                y_velocity=random.gauss(float(y[k+f]), float(y[k+f])*y_sd)
                x_pos=x_pos+x_velocity*(float(timestep)) #Random variation in wind at each time step
                y_pos=y_pos+y_velocity*float(timestep) #Random variation in wind at each time step
                if highest_windspeed[-1]<math.sqrt(math.pow(x_velocity, 2)+math.pow(y_velocity, 2)):
                    highest_windspeed[-1]=math.sqrt(math.pow(x_velocity, 2)+math.pow(y_velocity, 2))
            elapsed_time+=timestep
            prev_elev=elev[k+f]
            if k+f>len(t):
                break
            
        #Convert final position to lat/lon 
        lat2=lat1+(y_pos/111325) 
        lon2=lon1+(x_pos/(111325*(math.cos((math.pi/180)*(lat2+lat1)/2))))
        kml_file.write("<Placemark>\n")
        kml_file.write("<Name>Jake</Name>\n")
        kml_file.write("<visibility>1</visibility>\n")
        kml_file.write("<styleUrl>#msn_open-diamond</styleUrl>\n")
        kml_file.write("<Point>\n")
        kml_file.write("<altitudeMode>clampToGround</altitudeMode>\n")
        kml_file.write("<coordinates>"+str(lon2)+","+str(lat2)+",0</coordinates>\n")
        kml_file.write("</Point>\n")
        kml_file.write("</Placemark>\n")
        final_x.append(x_pos/1600)
        final_y.append(y_pos/1600)
    kml_file.write("</Document>\n")
    kml_file.write("</kml>")    
    kml_file.close()

    magnitude=[]
    for k in range(len(final_x)):
        magnitude.append(math.sqrt(final_x[k]*final_x[k]+final_y[k]*final_y[k]))
    print("Minimum flight distance: "+str(round(min(magnitude)))+ " miles")
    print("Maximum flight distance: "+str(round(max(magnitude)))+ " miles")
    print("Average flight distance: "+str(round(sum(magnitude)/len(magnitude)))+ " miles")
    print("*"*10)
    print("Lowest maximum windspeed: "+str(round(min(highest_windspeed)*3.2808399/5280*3600))+" mph")
    print("Highest maximum windspeed: "+str(round(max(highest_windspeed)*3.2808399/5280*3600))+" mph")
    print("Average maximum windspeed: "+str(round(sum(highest_windspeed)/len(highest_windspeed)*3.2808399/5280*3600))+" mph")
    print("*"*10)
    print("Maximum burst elevation:  "+str(round(max(max_elevs)*3.2808399))+" ft (maximum allowed by model is "+str(round(40000*3.2808399))+" ft")
    print("Minimum burst elevation:  "+str(round(min(max_elevs)*3.2808399))+" ft")
    print("Average burst elevation:  "+str(round(sum(max_elevs)/len(max_elevs)*3.2808399))+" ft")
    print(elapsed_time)
    return(0)


launch_point=[launch_lat, launch_lon]
forecast=0
#foo=monte_carlo_balloon_flight([launch_lat, launch_lon], flight_data, model_parameters)
model_format="list" #kml or list
launch_time=get_most_recent_model_time() #Get the most recent model off of the UW website
print("Forecast time (GMT): " +str(launch_time))
#foo, weather=get_trajectory_from_internet(launch_lat, launch_lon, 30000, launch_time, forecast, "kml")
foo, weather=get_trajectory_from_internet(launch_lat,launch_lon, 30000, launch_time, forecast, model_format)
#if model_format=="list":
#    foo=monte_carlo_balloon_flight(launch_point, weather, model_parameters)
#weather=read_uw_model(open("C:\\bovine_aerospace\\test_trajectory.txt", "r"))
ascent_model=build_model(weather, [0, 0, 0, 0])
descent_model=ascent_model
model_resolution=100
display=True
forward=True
mb=1
L=2.5
po=1.06
ma=0.02897
mp=0.004002602
mi=1
dcd=1.5
rp=0.3
rb=3.05
trials=1
file_name=str(launch_time)+"_point_"+str(round(launch_point[0], 3))+"_"+str(round(launch_point[1], 3))+".txt"
print("Get estimated landing point for descent model")
landing_point, xy=launch(launch_point, launch_time, ascent_model, descent_model, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, "initial_trajectory_"+file_name)
foo, descent_weather=get_trajectory_from_internet(landing_point[0],landing_point[1], 30000, launch_time, forecast, model_format)
descent_model=build_model(descent_weather, [0, 0, 0, 0])
landing_point, xy=launch(launch_point, launch_time, ascent_model, descent_model, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, file_name)
build_flight_kml(file_name)
#mc2(launch_point, launch_time, weather, model_resolution, forward, mb, L, po, ma, mp, mi, dcd, rp, rb, 1000)


