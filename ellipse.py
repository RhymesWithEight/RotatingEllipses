import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import unicodedata

# user generated values
print("Let's get our points.\n")
A = [float(item) for item in input("Enter the x and y values of the first point, with a comma separating them: ").split(',')]
B = [float(item) for item in input("Enter the x and y values of the second point, with a comma separating them: ").split(',')]
C = [float(item) for item in input("Enter the x and y values of the third point, with a comma separating them: ").split(',')]
I = int(input("How many ellipses iterations do you want? (Recommended 3-10): "))
d = float(input("As a decimal, how far do you want your curve to move along an ellipse each iteration? Must be between 0 and 2" + unicodedata.lookup("GREEK SMALL LETTER PI") + "(about 6.28). (Recommend less than 1):"))


N = 64 # resolution of ellipses
ms = 1 # speed for steps in animation, in seconds
dt = d / (2 * np.pi) # proportion d to entire circle

# some easy color references
b = 'b' # blue
r = 'r' # red
g = 'g' # green
k = 'k' # black
o = '#AAAAAA' # grey
l = '#EEEEEE' # light grey

def midpoint(H, K): # find the midpoint between two points
	return ((H[0]+K[0])/2, (H[1]+K[1])/2)

def distance(H, K): # find the distance between two points
	return np.sqrt((H[0]-K[0])**2+(H[1]-K[1])**2)

def theta_find(H, K): # find the angle between two points
	if H[0] == K[0]:
		return np.pi/2
	else:
		return np.arctan((H[1]-K[1])/(H[0]-K[0]))

def major(H, K, J): # length of the major axis for an ellipse with foci at H & K, containing the point J
	return (distance(H, J) + distance(K, J))/2

def ellipse(H, K, J, n): # create the x and y parameters of the ellipse, time n
	Oj = midpoint(H, K)
	Aj = major(H, K, J)
	Cj = distance(Oj, H)
	Bj = np.sqrt(Aj**2-Cj**2)
	Tj = theta_find(H, K)
	Xj = Aj*np.cos(Tj)*np.cos(n)-Bj*np.sin(Tj)*np.sin(n)+Oj[0]
	Yj = Aj*np.sin(Tj)*np.cos(n)+Bj*np.cos(Tj)*np.sin(n)+Oj[1]
	return [Xj, Yj]

def time_finder(H, K, J): # find the time J is generated with the parametric equations created with ellipse()
	Oj = midpoint(H, K)
	Aj = major(H, K, J)
	Cj = distance(Oj, H)
	Bj = np.sqrt(Aj**2-Cj**2)
	Tj = theta_find(H, K)
	
	# find cosine of t in order to identify it
	cos_t = 1/Aj*((J[0]-Oj[0])*np.cos(Tj)+(J[1]-Oj[1])*np.sin(Tj))
	sin_t = 1/Bj*(-(J[0]-Oj[0])*np.sin(Tj)+(J[1]-Oj[1])*np.cos(Tj))
	
	if sin_t == 0:
		return 0
	elif sin_t > 0:
		return np.arccos(cos_t)
	else:
		return -np.arccos(cos_t)

# prepare the plot
fig, ax = plt.subplots()

# set a bunch of variables needed to plot. x and y data lists for the ellipses, lines, points, etc.
Aplot, = ax.plot(A[0], A[1], marker=".", markersize=10, color=b)
Bplot, = ax.plot(B[0], B[1], marker=".", markersize=10, color=g)
Cplot, = ax.plot(C[0], C[1], marker=".", markersize=10, color=r)
A_graph, = ax.plot([],[], color=b)
B_graph, = ax.plot([],[], color=g)
C_graph, = ax.plot([],[], color=r)
A_line, = ax.plot([],[], color=k, linewidth=3)
B_line, = ax.plot([],[], color=k, linewidth=3)
C_line, = ax.plot([],[], color=k, linewidth=3)
A_new = A.copy()
B_new = B.copy()
C_new = C.copy()
A_moving = A.copy()
B_moving = B.copy()
C_moving = C.copy()
xdataC, ydataC = [], []
xdataA, ydataA = [], []
xdataB, ydataB = [], []
xdataLineC, ydataLineC = [], []
xdataLineA, ydataLineA = [], []
xdataLineB, ydataLineB = [], []

# set starting parameters
def init():
	t = np.linspace(0, 2*np.pi, N)
	
	max_xA = max(ellipse(B,C,A,t)[0])
	max_xB = max(ellipse(C,A,B,t)[0])
	max_xC = max(ellipse(A,B,C,t)[0])
	max_yA = max(ellipse(B,C,A,t)[1])
	max_yB = max(ellipse(C,A,B,t)[1])
	max_yC = max(ellipse(A,B,C,t)[1])
	
	min_xA = min(ellipse(B,C,A,t)[0])
	min_xB = min(ellipse(C,A,B,t)[0])
	min_xC = min(ellipse(A,B,C,t)[0])
	min_yA = min(ellipse(B,C,A,t)[1])
	min_yB = min(ellipse(C,A,B,t)[1])
	min_yC = min(ellipse(A,B,C,t)[1])
	
	max_x = max(max_xA, max_xB, max_xC)
	max_y = max(max_yA, max_yB, max_yC)
	min_x = min(min_xA, min_xB, min_xC)
	min_y = min(min_yA, min_yB, min_yC)
	
	ax.set_xlim(min_x, max_x, auto=True)
	ax.set_ylim(min_y, max_y, auto=True)

# draw all, in order
def draw(frame):
	drawA(frame)
	drawB(frame)
	drawC(frame)

# draw each ellipse
def drawA(frame):
	time_A_new = time_finder(B_new,C_new,A_new)
	frame_adjust = np.mod(frame, 2)
	if frame_adjust == 0 and frame != 0:
		A_new.clear()
		A_new.append(A_moving[0])
		A_new.append(A_moving[1])
		xdataDepA, ydataDepA = [], []
		xdataDepA.extend(xdataA)
		ydataDepA.extend(ydataA)
		A_dep, = ax.plot(xdataA, ydataA,color=l, linewidth=1, zorder=1)
		xdataA.clear()
		ydataA.clear()
		return A_dep
	elif frame_adjust < 1:
		A_graph.set_color(b)
		frame_adjust = 2*np.pi*frame_adjust
		xdataA.append(ellipse(B_new,C_new,A_new,frame_adjust+time_A_new-2*np.pi/N)[0])
		ydataA.append(ellipse(B_new,C_new,A_new,frame_adjust+time_A_new-2*np.pi/N)[1])
		A_graph.set_data(xdataA, ydataA)
		return A_graph,
	elif frame_adjust == 1:
		frame_adjust = 2*np.pi*frame_adjust
		xdataA.append(ellipse(B_new,C_new,A_new,frame_adjust+time_A_new-2*np.pi/N)[0])
		ydataA.append(ellipse(B_new,C_new,A_new,frame_adjust+time_A_new-2*np.pi/N)[1])
		xdataA.append(ellipse(B_new,C_new,A_new,time_A_new-2*np.pi/N)[0])
		ydataA.append(ellipse(B_new,C_new,A_new,time_A_new-2*np.pi/N)[1])
		xdataA.append(ellipse(B_new,C_new,A_new,time_A_new)[0])
		ydataA.append(ellipse(B_new,C_new,A_new,time_A_new)[1])
		A_graph.set_data(xdataA, ydataA)
		return A_graph,
	else:
		frame_adjust = 2*np.pi*frame_adjust
		A_graph.set_color(o)
		t_current = time_finder(B_new,C_new,A_new) + (frame_adjust - 2*np.pi)*dt
		A_moving.clear()
		A_moving.append(ellipse(B_new,C_new,A_new,t_current)[0])
		A_moving.append(ellipse(B_new,C_new,A_new,t_current)[1])
		xdataLineA.append(A_moving[0])
		ydataLineA.append(A_moving[1])
		A_line.set_data(xdataLineA, ydataLineA)
		Aplot.set_data(A_moving[0],A_moving[1])
		return A_line,

def drawB(frame):
	time_B_new = time_finder(C_new,A_new,B_new)
	frame_adjust = np.mod(frame, 2)
	if frame_adjust == 0 and frame != 0:
		B_new.clear()
		B_new.append(B_moving[0])
		B_new.append(B_moving[1])
		xdataDepB, ydataDepB = [], []
		xdataDepB.extend(xdataB)
		ydataDepB.extend(ydataB)
		B_dep, = ax.plot(xdataB, ydataB,color=l, linewidth=1, zorder=1)
		xdataB.clear()
		ydataB.clear()
		return B_dep,
	elif frame_adjust < 1:
		B_graph.set_color(g)
		frame_adjust = 2*np.pi*frame_adjust
		xdataB.append(ellipse(C_new,A_new,B_new,frame_adjust+time_B_new-2*np.pi/N)[0])
		ydataB.append(ellipse(C_new,A_new,B_new,frame_adjust+time_B_new-2*np.pi/N)[1])
		B_graph.set_data(xdataB, ydataB)
		return B_graph,
	elif frame_adjust == 1:
		frame_adjust = 2*np.pi*frame_adjust
		xdataB.append(ellipse(C_new,A_new,B_new,frame_adjust+time_B_new-2*np.pi/N)[0])
		ydataB.append(ellipse(C_new,A_new,B_new,frame_adjust+time_B_new-2*np.pi/N)[1])
		xdataB.append(ellipse(C_new,A_new,B_new,time_B_new-2*np.pi/N)[0])
		ydataB.append(ellipse(C_new,A_new,B_new,time_B_new-2*np.pi/N)[1])
		xdataB.append(ellipse(C_new,A_new,B_new,time_B_new)[0])
		ydataB.append(ellipse(C_new,A_new,B_new,time_B_new)[1])
		B_graph.set_data(xdataB, ydataB)
		return B_graph,
	else:
		frame_adjust = 2*np.pi*frame_adjust
		B_graph.set_color(o)
		t_current = time_finder(C_new,A_new,B_new) + (frame_adjust - 2*np.pi)*dt
		B_moving.clear()
		B_moving.append(ellipse(C_new,A_new,B_new,t_current)[0])
		B_moving.append(ellipse(C_new,A_new,B_new,t_current)[1])
		xdataLineB.append(B_moving[0])
		ydataLineB.append(B_moving[1])
		B_line.set_data(xdataLineB, ydataLineB)
		Bplot.set_data(B_moving[0],B_moving[1])
		return B_line,

def drawC(frame):
	time_C_new = time_finder(A_new,B_new,C_new)
	frame_adjust = np.mod(frame, 2)
	if frame_adjust == 0 and frame != 0:
		C_new.clear()
		C_new.append(C_moving[0])
		C_new.append(C_moving[1])
		xdataDepC, ydataDepC = [], []
		xdataDepC.extend(xdataC)
		ydataDepC.extend(ydataC)
		C_dep, = ax.plot(xdataC, ydataC,color=l, linewidth=1, zorder=1)
		xdataC.clear()
		ydataC.clear()
		return C_dep,
	elif frame_adjust < 1:
		C_graph.set_color(r)
		frame_adjust = 2*np.pi*frame_adjust
		xdataC.append(ellipse(A_new,B_new,C_new,frame_adjust+time_C_new-2*np.pi/N)[0])
		ydataC.append(ellipse(A_new,B_new,C_new,frame_adjust+time_C_new-2*np.pi/N)[1])
		C_graph.set_data(xdataC, ydataC)
		return C_graph,
	elif frame_adjust == 1.0:
		frame_adjust = 2*np.pi*frame_adjust
		xdataC.append(ellipse(A_new,B_new,C_new,frame_adjust+time_C_new-2*np.pi/N)[0])
		ydataC.append(ellipse(A_new,B_new,C_new,frame_adjust+time_C_new-2*np.pi/N)[1])
		xdataC.append(ellipse(A_new,B_new,C_new,time_C_new-2*np.pi/N)[0])
		ydataC.append(ellipse(A_new,B_new,C_new,time_C_new-2*np.pi/N)[1])
		xdataC.append(ellipse(A_new,B_new,C_new,time_C_new)[0])
		ydataC.append(ellipse(A_new,B_new,C_new,time_C_new)[1])
		C_graph.set_data(xdataC, ydataC)
		ax.relim()
		ax.autoscale_view()
		return C_graph,
	else:
		frame_adjust = 2*np.pi*frame_adjust
		C_graph.set_color(o)
		t_current = time_finder(A_new,B_new,C_new) + (frame_adjust - 2*np.pi)*dt
		C_moving.clear()
		C_moving.append(ellipse(A_new,B_new,C_new,t_current)[0])
		C_moving.append(ellipse(A_new,B_new,C_new,t_current)[1])
		xdataLineC.append(C_moving[0])
		ydataLineC.append(C_moving[1])
		C_line.set_data(xdataLineC, ydataLineC)
		Cplot.set_data(C_moving[0],C_moving[1])
		return C_line,

# actually animate
draw_ellipses = animation.FuncAnimation(fig, draw, frames=np.linspace(0, 2*I, 2*I*N+1), init_func=init, interval=ms, repeat=False)

# show the result
plt.show()