# BSON encoding

BSON_enoder.py is a Python script for delta encoding GeoJSON files and converting them to BSON - the binary file format of MongoDB. By converting the files to BSON as well as applying delta encoding and using GZIP, a file size reduction of almost 90% can be achieved, while retaining full coordinate precision. 
