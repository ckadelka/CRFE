# -*- coding: utf-8-sig -*-

#Only recognizes miRNAs that contain mir in the name

OUTPUT_INTO_FILE=True
KNOCKOUT_MIRNAS=False

filename='updaterules61113mir1792minus'

add_info="DNA Mismatch Repair Pathway. REU Summer 2013 at VBI"

f = open(filename+'.txt','r')
text=f.read()

text=text.replace(" ","")

text=text.replace("[","(")
text=text.replace("] ",")")
text=text.replace("{","(")
text=text.replace("} ",")")

ADAM=0
MATLAB=1
PYTHON=0

if ADAM:
    sign=["", ""]
elif MATLAB:
    sign=["(",")"]
elif PYTHON:
    sign=["[","]"]

tvec=text.splitlines()

#Remove empty lines
while True:
    try:
        tvec.remove('')
    except:
        break
        
n=len(tvec)
var=["" for i in xrange(n)]
s=["" for i in xrange(n)]
outvec=["" for i in xrange(n)]
for i in xrange(n):
    var[i]=tvec[i][0:tvec[i].find("=")]

for i in xrange(n,0,-1):
    ind=tvec[i-1].find('=')
    ind2=tvec[i-1].find('=',ind+1)
    s[i-1]=tvec[i-1][ind+1:ind2]
    tvec[i-1]=tvec[i-1][ind2+1:]
    #print "$"+var[i-1]+"$","$"+tvec[i-1]+"$"
    if 'mir' in var[i-1] and KNOCKOUT_MIRNAS==True:
        tvec[i-1]=var[i-1]
        s[i-1]='1'
    for j in xrange(n,0,-1):
        tvec[i-1]=tvec[i-1].replace(var[j-1],"x"+sign[0]+str(j)+sign[1])
    
#for i in xrange(n):
#    tvec[i]=tvec[i].replace("!","~")
#    tvec[i]=tvec[i].replace("OR","|")
#    tvec[i]=tvec[i].replace("AND","&")


invar=[0 for i in xrange(n)]
varF=[]
for i in xrange(n):      
    count=0
    for j in xrange(n):
        if tvec[i].find("x"+sign[0]+str(j+1)+sign[1]) > -1:
            tvec[i]=tvec[i].replace("x"+sign[0]+str(j+1)+sign[1],"x"+sign[0]+str(count+1)+sign[1])
            count+=1
            if count>len(varF):
                varF.append([-1 for k in xrange(n)])
            varF[count-1][i]=j+1
    invar[i]=count

varFstr=""
for i in xrange(len(varF)):
    helper=str(varF[i])
    helper=helper[1:-1]
    varFstr+=helper+";"
varFstr=varFstr[:-1]

s_int=[]
for i in xrange(len(s)):
    s_int.append(int(s[i]))
    

ai=add_info.splitlines()

extension="_miRNAs_all_0" if KNOCKOUT_MIRNAS==True else ""

output="function [res,text,varF,s,names] = "+filename+extension+"(x,i,dummy)\n"+"\n" 
output+="text ='"+add_info+"';\n"+"\n" 
output+="varF = ["+varFstr+"];\n"+"\n"     
output+="s = "+str(s_int)+";\n"+"\n"     
helper=str(var)
output+="names = {"+helper[1:-1]+"};\n"+"\n"
output+="if nargin==3"+"\n"
output+="\tres = "+str(invar)+";"+"\n"
output+="else"+"\n"
for i in xrange(n):
    outvec[i] ="res = "+tvec[i]
    if MATLAB:
        outvec[i]=outvec[i]+";"
    if i==0:
        output+="\tif i=="+str(i+1)+"\n"
    else:
        output+="\telseif i=="+str(i+1)+"\n"
    output+="\t\t"+outvec[i]+"\n"
output+="\tend"+"\n"
output+="end"+"\n"

if OUTPUT_INTO_FILE==True:
    g=open(filename+extension+'.m','w')
    g.write(output)
    g.close()
else:
    print output