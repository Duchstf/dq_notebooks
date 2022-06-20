import uproot
import awkward as ak
import matplotlib.pyplot as plt

def vertex_std(root_file, plot=False):
	# Read the file
	f = uproot.open(root_file) #'/home/caspian/data/Brem_2.750000_z500_600_eps_-6.8.root'
	x = f['Events'].arrays()
	EventID
	if len(x) < 500:  # the results are not useful if we don't have a lot of events
		raise Exception

	# Filter hit_pos to only include data with DP Hodo hits
	hit_y = x.hit_pos[ (55 <= x.hit_detID) & (x.hit_detID <= 62) ]
	filtered_detID = x.hit_detID[ (55 <= x.hit_detID) & (x.hit_detID <= 62) ]

	# Filter to require 4 hits (used to draw two lines)
	have_4_hits_event_mask = ak.sum(hit_y!=None,axis=1) == 4
	hit_y = hit_y[ have_4_hits_event_mask ]
	filtered_detID = filtered_detID[ have_4_hits_event_mask ]

	# Filter to require hits to stay in same quadrant
	same_quadrant_mask = (filtered_detID[:,0] + 4 == filtered_detID[:,2]) & (filtered_detID[:,1] + 4 == filtered_detID[:,3])
	# Filter to require hits to stay in same quadrant, but be in same or different quadrant
	#same_quadrant_mask = (filtered_detID[:,0] == filtered_detID[:,1]) & (filtered_detID[:,0] + 4 == filtered_detID[:,2]) & (filtered_detID[:,1] + 4 == filtered_detID[:,3])
	print('Acceptance :',ak.sum(same_quadrant_mask)/len(same_quadrant_mask), root_file)
	# Uncomment the line below to make the acceptance plot (rather than std error plot)
	#return ak.sum(same_quadrant_mask)/len(same_quadrant_mask)
	hit_y = hit_y[ same_quadrant_mask ]
	filtered_detID = filtered_detID[ same_quadrant_mask ]

	# Calculate the slopes of the lines
	m1=(hit_y[:,2]-hit_y[:,0])/(1471.0-798.0)
	m2=(hit_y[:,3]-hit_y[:,1])/(1471.0-798.0)

	# Calculate the z coordinate of the intersection
	x1f,y1=798,hit_y[:,0]
	x2f,y2=798,hit_y[:,1]
	vertex_from_track = (-m2*x2f+y2-y1+m1*x1f ) / (m1-m2)

	truth_vertex = x.truthtrack_z_vtx[:,0]
	truth_vertex = truth_vertex[have_4_hits_event_mask] # apply the same mask on events
	truth_vertex = truth_vertex[same_quadrant_mask] # apply the same mask on events

	# Plot
	error = vertex_from_track - truth_vertex
	error = error[error > -50]
	error = error[error < 50]
	import numpy as np
	import scipy.stats
	#print('std',np.std(error), root_file)
	if plot:
		plt.plot(np.linspace(-50,50),scipy.stats.norm.pdf(np.linspace(-50,50),0,np.std(error)),label='Gaussian fit\nstd='+round(np.std(error),2)+' cm')
		plt.legend()
		plt.xlabel('Distribution of error in z vertex coordinate (cm)')
		plt.hist(error,bins=30,density=True)
		plt.title('Error in Calculated Vertex z Position')
		plt.savefig('output.png')
	else:
		return np.std(error)





import glob

eps=[]  # coupling constant
mass=[]  
mech = "Brem"  # A' generate mechanisim bremsstrahlung
std=[]

for root_file in glob.glob('/home/caspian/data/{}*.root'.format(mech)):
	info      = root_file.split('_')
	mass_info = float(info[1])
	eps_info  = float(info[-1][:-5])
	try:
		std.append(vertex_std(root_file,plot=False))
		mass.append(mass_info)
		eps.append(eps_info)
	except Exception:
		print('failed',root_file)

# following https://stackoverflow.com/questions/25581396/2d-plot-with-matplotlib
plt.scatter(mass, eps, c=std ,s=200, marker='s',cmap='Spectral_r',linewidths=0)
plt.colorbar()
plt.clim(0,30)
#plt.clim(0,1)
plt.title('Calculated vertex z coordinate error for {} with requiring both roads to be in the same quadrant'.format(mech))
#plt.title(f'Acceptance for {mech} with no requirements on roads')
#plt.title(f'Acceptance for {mech} with requiring roads to be in different quadrants')
#plt.title(f'Acceptance for {mech} with requiring both roads to be in the same quadrant')
plt.ylabel('Coupling Strength log(eps)')
plt.xlabel('Mass (GeV)')
plt.show()





