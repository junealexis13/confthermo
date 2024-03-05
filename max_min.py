nwindow = int(input("write nwindow:"))
f1=open("max_min_apo.dat","r")
f2=open("max_min_holo.dat","r")
f3=open("max_min.dat","w")
for i in range(nwindow):
    x,y = f1.readline().split()
    x=float(x)
    y=float(y)
    x1,y1 = f2.readline().split()
    x1=float(x1)
    y1=float(y1)
    if (x1>x) and (y1<y):
        print ("%12.8f  %12.8f"%(x1,y1),file=f3)
        print ("%12.8f  %12.8f"%(x1,y1),file=f3)
    elif (x1<x) and (y1<y):
        print ("%12.8f  %12.8f"%(x,y1),file=f3)
        print ("%12.8f  %12.8f"%(x,y1),file=f3)
    elif (x1>x) and (y1>y):
        print ("%12.8f  %12.8f"%(x1,y),file=f3)
        print ("%12.8f  %12.8f"%(x1,y),file=f3)
    else:
        print ("%12.8f  %12.8f"%(x,y),file=f3)
        print ("%12.8f  %12.8f"%(x,y),file=f3)
f1.close()
f2.close()
f3.close()

f4=open("max_min.dat","r")
f5=open("max_min_convert.dat","w")
for i in range(2*nwindow):
    a,b = f4.readline().split()
    a=float(a)
    b=float(b)
    print(a,file=f5)
    print(b,file=f5)
f4.close()
f5.close()
