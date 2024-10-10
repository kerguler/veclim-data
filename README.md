# VEClim's data server

# Example usage
 Run the server with the following:

 -   `python3 veclim-python-server.py &> log &`

 Try the following requests:

 - To get the decadal mean on a single day:
    
    `http://localhost:${VEC_PORTE}/?lon=33.3&lat=35.0&date=2020-02-28`

 - To get the decadal mean over a period:

    `http://localhost:${VEC_PORTE}/?lon=33.3&lat=35.0&dates=2020-02-28:2020-03-07`

 - To get the decadal mean for each day over a period:

    `http://localhost:${VEC_PORTE}/?lon=33.3&lat=35.0&dates=2020-02-28:2020-03-07&opr=ts`

Please note that you will need to configure a .env file in the docker/VEClim/veclim-server directory as follows:
 - DIR_DATA=The internal data directory containing the environmental datasets and simulation outputs
 - DIR_PYTHON=The dist-paackages directory of the python installation inside the docker image
 - VEC_DATA=The external data directory containing the environmental datasets and simulation outputs
 - VEC_HOST=Desired host IP of the python server
 - VEC_PORT=Desired internal port of the python server
 - VEC_PORTE=Desired external port of the python server