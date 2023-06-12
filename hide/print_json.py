labels=["Topping","Slices"]
types=["string","number"]
keys=["Mushrooms","Salami","Ham"]
values=[3,1,2]

#jsonDATA = "{\n\"cols\": [\n"
#for (label,typ) in zip(labels,types):
#    jsonDATA += "{\"label\":\""+label+"\",\"type\":\""+typ+"\"},\n"
#jsonDATA=jsonDATA[:-2]+"\n],\n\"rows\": ["
#for (key,value) in zip(keys,values):
#    jsonDATA += "{\"c\":[{\"v\":\""+key+"\",\"f\":null},{\"v\":"+str(value)+",\"f\":null}]},\n"
#jsonDATA=jsonDATA[:-2]+"\n]\n}"

jsonDATA = "{\"cols\": ["
for (label,typ) in zip(labels,types):
    jsonDATA += "{\"label\":\""+label+"\",\"type\":\""+typ+"\"},"
jsonDATA=jsonDATA[:-1]+"],\"rows\": ["
for (key,value) in zip(keys,values):
    jsonDATA += "{\"c\":[{\"v\":\""+key+"\",\"f\":null},{\"v\":"+str(value)+",\"f\":null}]},"
jsonDATA=jsonDATA[:-1]+"]}"

print jsonDATA