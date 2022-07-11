







function etr = estimated_time_remaining(startpoint,endpoint,currentpoint,timetakensofar)



totalrange = (endpoint-startpoint+1);

avgtime = (currentpoint-startpoint+1)/timetakensofar;

esttotaltime = avgtime * totalrange;


etr = esttotaltime - timetakensofar;

etr = etr / 60; % time remaining in minutes


























