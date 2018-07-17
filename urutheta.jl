# Read file.
matrix=readdlm("temp_axisym.txt")
# Allocate cartesian speeds: u-speed
u=matrix[:,3]
v=matrix[:,4]
r=matrix[:,1]
theta=matrix[:,2]

H=matrix[:,6]
siz=size(v)
v_radial=ones(siz[1])
v_azi=ones(siz[1])
xi=ones(siz[1])
yi=ones(siz[1])
for i in range(1,siz[1])
v_radial[i]=round(u[i]*cos(theta[i])+v[i]*sin(theta[i]),3)
v_azi[i]=round(-u[i]*sin(theta[i])+v[i]*cos(theta[i]),3)
xi[i]=round(r[i]*cos(theta[i]),3)
yi[i]=round(r[i]*sin(theta[i]),3)
end
print(size(matrix),size(v_azi))
matrix = cat(1,matrix', v_azi')
matrix=cat(1,matrix,v_radial')
matrix = cat(1,matrix, xi')
matrix = cat(1,matrix, yi')
writedlm("tempjulia.txt",matrix')
