# August 11, 2023 Robert Holder
# This code estimates the total runtime for a laser ablation experiment on the NWR193 excimer laser ablation system.
# To use this script:
# 1) paste your laser experiment file (.LAX) into the appropriate working folder
# 2) enter the values for the washout delay and warmup time as well as the laser experiment file name into the input file
# 3) run the script however you run python on your computer
################################################################################

# import libraries
import re #regular expressions for text searching
import numpy as np #numpy to use numpy arrays, which are a bit easier than native lists


def estimate_runtime(filename, warmup, washout):
    stagespeed = 1300 # um/s estimated from watching laser move during a run. should be decently accurate
    waiting = warmup + washout

    # Try different encodings.
    encodings = ['utf-8', 'latin-1', 'iso-8859-1', 'utf-16']
    
    for encoding in encodings:
        try:
            with open(filename, 'r', encoding=encoding) as f:
                text = f.read()
            break  # If successful, exit the loop
        except UnicodeDecodeError:
            continue  # If unsuccessful, try the next encoding
    else:
        # If all encodings fail, raise an error
        raise ValueError("Unable to decode the file with any of the attempted encodings")

    # spots or line scans?
    spotsorlines = re.search('<Type>(.+?)</Type>',text)

    if spotsorlines[1] == 'Spot':
        x = re.findall('<PointList>[\S\s]*?<X>([0-9]+.[0-9]+)</X>',text)
        y = re.findall('<PointList>[\S\s]*?<Y>([0-9]+.[0-9]+)</Y>',text)
        Xspot = np.array([float(i) for i in x])
        Yspot = np.array([float(i) for i in y])

        # total time for spot ablations
        spotdwell1 =  re.search('<DwellTime>(.+?)</DwellTime>',text)
        spotdwell = float(spotdwell1[1])
        spotduration = round(len(Xspot)*spotdwell / 60 / 60,2)

        # estimate the total duration associated with warmup and washout delay
        preablation = re.search('<EnablePreAblationPass>(.+?)</EnablePreAblationPass>',text)
        if preablation[1] == "True":
            spotwaiting = 2* round(waiting*len(Xspot) / 60 / 60,2) # unit = hrs
        else:
            spotwaiting = round(waiting*len(Xspot) / 60 / 60,2) # unit = hrs

        # calculate the total time required to move the stage from one spot to the next
        spotspacing = []
        for i in range(len(Xspot)-1):
            spotspacing.append( abs(Xspot[i+1]-Xspot[i]) + abs(Yspot[i+1]-Yspot[i]) )
        spotspacing = np.array(spotspacing)
        spotmove = round(sum(spotspacing / stagespeed) / 60 / 60,2) # unit = hrs

        # calculate total runtime by summing the time estimates above.
        spotruntime = round(spotduration+spotmove+spotwaiting,2)
        spotruntimefudge = round(spotruntime*1.10,2)

        print(  'total number of spots           ',len(Xspot),
                '\n\ntotal analysis time              ',int(spotduration//1),'hrs',int(spotduration%1*60),'min',
                '\ntotal movement time              ',int(spotmove//1),'hrs',int(spotmove%1*60),'min',
                '\ntotal warmup+washout             ',int(spotwaiting//1),'hrs',int(spotwaiting%1*60),'min',
                '\ntotal estimated runtime          ',int(spotruntime//1),'hrs',int(spotruntime%1*60),'min',
                '\nwith an extra 10% fudge factor   ',int(spotruntimefudge//1),'hrs',int(spotruntimefudge%1*60),'min',)
        input("\n\n\npress enter to close")


    else:
        # find the x and y coordinates at the start and end of each line and write them to arrays
        x1 = re.findall('<PointList>[\S\s]*?<X>([0-9]+.[0-9]+)</X>',text)
        x2 = re.findall('<PointList>[\S\s]*?<X>[\S\s]*?<X>([0-9]+.[0-9]+)</X>',text)
        Xstart = np.array([float(i) for i in x1])
        Xend = np.array([float(i) for i in x2])
        y1 = re.findall('<PointList>[\S\s]*?<Y>([0-9]+.[0-9]+)</Y>',text)
        y2 = re.findall('<PointList>[\S\s]*?<Y>[\S\s]*?<Y>([0-9]+.[0-9]+)</Y>',text)
        Ystart = np.array([float(i) for i in y1])
        Yend = np.array([float(i) for i in y2])

        # calculate the length of each scan (distance from start to end)
        linelength = ( (Xstart-Xend)**2 + (Ystart-Yend)**2 )**0.5

        # calculate the distance from the end of one scan to the start of the next
        linespacing = []
        for i in range(len(linelength)-1):
            linespacing.append( abs(Xstart[i+1]-Xend[i]) + abs(Ystart[i+1]-Yend[i]) )
        linespacing = np.array(linespacing)

        # extract the scan speed from the laser experiment file. Assumes the first scan has the same scanspeed as every other scan
        scanspeed = float(re.search('<BlastingProperties>[\S\s]*?<ScanSpeed>([0-9]+)</ScanSpeed>',text)[1])
        # calculate the total time required to complete all scans from scanspeed and scanlengths
        lineduration = round(sum(linelength / scanspeed) / 60 / 60,2) # unit = hrs

        # calculate the total time required to move the stage from the end of each scan to the start of the next
        linemove = round(sum(linespacing / stagespeed) / 60 / 60,2) # unit = hrs

        # estimate the total duration associated with warmup and washout delay
        linewaiting = round(waiting*len(Xstart) / 60 / 60,2) # unit = hrs

        # calculate total runtime by summing the time estimates above.
        lineruntime = round(lineduration+linemove+linewaiting,2)
        lineruntimefudge = round(lineruntime*1.1,2)

        print(  'total number of scans            ',len(Xstart),
                '\n\ntotal analysis time              ',int(lineduration//1),'hrs',int(lineduration%1*60),'min',
                '\ntotal movement time              ',int(linemove//1),'hrs',int(linemove%1*60),'min',
                '\ntotal warmup+washout             ',int(linewaiting//1),'hrs',int(linewaiting%1*60),'min',
                '\ntotal estimated runtime          ',int(lineruntime//1),'hrs',int(lineruntime%1*60),'min',
                '\nwith an extra 10% fudge factor   ',int(lineruntimefudge//1),'hrs',int(lineruntimefudge%1*60),'min',)
        input("\n\n\npress enter to close")




