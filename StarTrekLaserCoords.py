import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np

PI = 3.141592654
SECONDS_PER_POINT = 0.001
FPS = 15

def perform_interpolation(input_times, input_coord, output_coord):
    output_coord.clear()
    
    current_time = input_times[0];
    while current_time <= input_times[-1]:
        # find p0, p1, m0, and m1

        # find the index of the previous point
        k_0 = 0
        equal_done = False
        
        interpolated_point = 0

        for i in range(len(input_times)):
            if (input_times[i] == current_time):
                # just put the point in and be done with it
                interpolated_point = input_coord[i]
                equal_done = True;
                break
            if (input_times[i] > current_time):
                k_0 = i - 1
                break
        if (not equal_done):
            # find the index of the next point
            k_1 = k_0 + 1

            # find the slope of the previous point
            m_0 = 0
            if (k_0 == 0):
                # edge case when k_0 == 0
                m_0 = (input_coord[k_1] - input_coord[k_0]) / (input_times[k_1] - input_times[k_0])
            else:
                # k_0>0, so k_0-1 exists
                m_0 = 0.5 * ((input_coord[k_1] - input_coord[k_0]) / (input_times[k_1] - input_times[k_0]) + (input_coord[k_0] - input_coord[k_0 - 1]) / (input_times[k_0] - input_times[k_0 - 1]))

            # find the slope of the next point
            m_1 = 0
            if (k_1 == len(input_times) - 1):
                # edge case when k_1 == len(input_times) - 1
                m_1 = (input_coord[k_1] - input_coord[k_0]) / (input_times[k_1] - input_times[k_0])
            else:
                # k_1< len(input_times) - 1, so k_1+1 exists
                m_1 = 0.5 * ((input_coord[k_1 + 1] - input_coord[k_1]) / (input_times[k_1 + 1] - input_times[k_1]) + (input_coord[k_1] - input_coord[k_0]) / (input_times[k_1] - input_times[k_0]))

            # find the distance we need to interpolate as a double between 0 and 1, 0 being point 0, and 1 being point 1
            t = (current_time - input_times[k_0]) / (input_times[k_1] - input_times[k_0])

            # use t to get h00, h10, h01, and h11
            h00 = 2 * t * t * t - 3 * t * t + 1
            h10 = t * t * t - 2 * t * t + t
            h01 = -2 * t * t * t + 3 * t * t
            h11 = t * t * t - t * t

            # try renormalizing m_0 and m_1
            m_0 *= (input_times[k_1] - input_times[k_0])
            m_1 *= (input_times[k_1] - input_times[k_0])

            interpolated_point = h00 * input_coord[k_0] + h10 * m_0 + h01 * input_coord[k_1] + h11 * m_1

        # put in the interpolated coordinate
        output_coord.append(interpolated_point)

        # continue with the next time step
        current_time += SECONDS_PER_POINT
    return output_coord
        
def plot_path(verts = []):
    verts = [
       (0., 0.),   # P0
       (0.2, 1.),  # P1
       (1., 0.8),  # P2
       (0.8, 0.),  # P3
    ]
    
    colors = [
       (255, 255, 0.),   # P0
       (0., 255, 255),  # P1
       (255, 0., 255),  # P2
       (0., 0., 0.),  # P3
    ]
    
    xs, ys = zip(*verts)
    Xs, Ys = [], []
    
    rs, gs, bs = zip(*colors)
    Rs, Gs, Bs = [], [], []
    
    perform_interpolation([0,1,2,3], xs, Xs)
    perform_interpolation([0,1,2,3], ys, Ys)
    perform_interpolation([0,1,2,3], rs, Rs)
    perform_interpolation([0,1,2,3], gs, Gs)
    perform_interpolation([0,1,2,3], bs, Bs)
    
    
    Xs = np.clip(Xs, 0., 255.)
    Ys = np.clip(Ys, 0., 255.)
    Rs = np.clip(Rs, 0., 255.)
    Gs = np.clip(Gs, 0., 255.)
    Bs = np.clip(Bs, 0., 255.)
    
    Colors = np.array(list(zip(Rs, Gs, Bs)))
    plt.scatter(Xs, Ys, s=1, c=Colors/255)
    
    # write to file
    f = open("LaserCoordsOutput.txt", "w")
    
    f.write("Xs\n")
    for x in Xs:
        f.write(str(round(x)))
        f.write(",")
    f.write("\nYs\n")
    for y in Ys:
        f.write(str(round(y)))
        f.write(",")
    f.write("\nRs\n")
    for r in Rs:
        f.write(str(round(r)))
        f.write(",")
    f.write("\nGs\n")
    for g in Gs:
        f.write(str(round(g)))
        f.write(",")
    f.write("\nBs\n")
    for b in Bs:
        f.write(str(round(b)))
        f.write(",")
    
    plt.show()


def main():
    print("Hello World!")
    plot_path()
    
if __name__ == "__main__":
    main()