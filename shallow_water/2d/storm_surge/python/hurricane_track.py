#!/usr/bin/env python

import sys
import numpy as np

def date_to_time(date,start_day=0,start_hour=0):
    month = date.split('/')[0]
    day = date.split('/')[1]
    hour = date.split('/')[2][:-1] # Strips of Z at the end
    
    return ((int(day) - start_day) * 24 + (int(hour) - start_hour)) * 60**2
    

def convert_hurricane_track_data(data_path,out_path='hurricane_track.data',plot=False):
    # Data arrays
    latitude = []
    longitude = []
    time = []
    wind = []
    pressure = []

    # Load in hurricane track information
    hurricane_track_file = open(data_path,'r')
    
    # Header
    date = hurricane_track_file.readline()
    hurricane_name = hurricane_track_file.readline()
    hurricane_track_file.readline() 
    
    # Read in first line for time delta
    line = hurricane_track_file.readline()
    data_line = line.split()
    latitude.append(data_line[1])
    longitude.append(data_line[2])
    start_day = int(data_line[3].split('/')[1])
    start_hour = int(data_line[3].split('/')[2][:-1])
    time.append(date_to_time(data_line[3],start_day,start_hour))
    wind.append(data_line[4])
    pressure.append(data_line[5])

    for line in hurricane_track_file:
        data_line = line.split()
        latitude.append(data_line[1])
        longitude.append(data_line[2])
        time.append(date_to_time(data_line[3],start_day,start_hour))
        wind.append(data_line[4])
        pressure.append(data_line[5])

    hurricane_track_file.close()

    # Write out new file
    out_file = open(out_path,'w')
    
    # Write out header
    out_file.write("%s =: Number of time points\n" % len(time))
    out_file.write("\n")

    # Write out data in (time,lat,long,wind,pressure)
    for (i,t) in enumerate(time):
        out_file.write("%s %s %s %s %s\n" % (t,latitude[i],longitude[i],wind[i],pressure[i]))
        
    out_file.close()
    
    if plot:
        import matplotlib.pyplot as plt
        
        plt.figure(1)
        plt.plot(longitude,latitude,'o')
    
        fig = plt.figure(2)
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        ax1.plot(time,wind,'k')
        ax2.plot(time,pressure,'r')
    
        plt.show()
        
if __name__ == "__main__":
    if len(sys.argv) > 1:
        convert_hurricane_track_data(*sys.argv[1:])
    else:
        convert_hurricane_track_data('./track.data')