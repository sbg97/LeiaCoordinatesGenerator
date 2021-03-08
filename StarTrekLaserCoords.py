from PointInterpolatorAndOutputter import interpolateAndOutput
import math

SECONDS_PER_FRAME = 1./15.

PI = 3.141592654
explosioncoords = [(0,0), (1,4), (-5,8), (-2,1), (-8,0), (-3,-2), (-7,-7), (-1,-3)]
explosion2coords = [(0,0), (1,4), (-5,8), (-2,1), (-8,0), (-3,-2), (-7,-7), (-1,-3), (-1,-8), (2,-4)]

def shootDisruptor(verts, times, colors, start, end, percentage = 0.4, duration = 0.249):
    start_time = times[-1]
    t = start_time
    
    #define 2 points per frame
    while(t < start_time + duration):
        #insert first point
        t += SECONDS_PER_FRAME / 2.0
        times.append(t)
        
            #interpolate between start and end
        interp = max((1.0 + percentage) / duration * (t - start_time) - percentage, 0.0)
        verts.append(((1-interp)*start[0] + (interp)*end[0], (1-interp)*start[1] + (interp)*end[1]))
        
        colors.append((0,255,0))
        
        #insert second point
        t += SECONDS_PER_FRAME / 2.0
        times.append(t)
        
            #interpolate between start and end
        interp = min((1.0 + percentage) / duration * (t - start_time), 1.0)
        verts.append(((1-interp)*start[0] + (interp)*end[0], (1-interp)*start[1] + (interp)*end[1]))
        
        colors.append((0,255,0))
        
def shootPhaser(verts, times, colors, start, end1, end2, duration = 1.5):
    start_time = times[-1]
    t = start_time
    
    #define 2 points per frame
    while(t < start_time + duration):
        #interpolate between end1 and end2
        interp = (t - start_time) / duration
        
        #insert first point
        t += SECONDS_PER_FRAME / 2.0
        times.append(t)
        
        verts.append(start)
        
        colors.append((255,0,0))
        
        #insert second point
        t += SECONDS_PER_FRAME / 2.0
        times.append(t)
        
        verts.append(((1-interp)*end1[0] + (interp)*end2[0], (1-interp)*end1[1] + (interp)*end2[1]))
        
        colors.append((255,0,0))
        
def littleExplosion(verts, times, colors, start, duration = 0.1295, rotation = 100.0 / 180.0 * PI, scale = 4.0):
    start_time = times[-1]
    t = start_time
    
    # define N points per frame where N is len(explosioncoords)
    while(t < start_time + duration):
        for v in explosioncoords:
            t += SECONDS_PER_FRAME / len(explosioncoords)
            times.append(t)
            
            # 2D translate, rotate, and scale the points of v
            verts.append((start[0] + scale * (v[0]*math.cos(rotation) + v[1]*math.sin(rotation)), (start[1] + scale * (v[1]*math.cos(rotation) - v[0]*math.sin(rotation)))))
            
            colors.append((255,0,0))
            
def bigExplosion(verts, times, colors, start, duration = 1.484, rotation = -90.0 / 180.0 * PI, scale1 = 2.0, scale2 = 5.0):
    start_time = times[-1]
    t = start_time
    
    # define N points per frame where N is len(explosioncoords)
    while(t < start_time + duration):
        for v in explosion2coords:
            
            #interpolate between scale1 and scale2
            interp = (t - start_time) / duration
            scale = scale1 * (1.0 - interp) + scale2 * interp
            
            t += SECONDS_PER_FRAME / len(explosion2coords)
            times.append(t)
            
            # 2D translate, rotate, and scale the points of v
            verts.append((start[0] + scale * (v[0]*math.cos(rotation) + v[1]*math.sin(rotation)), (start[1] + scale * (v[1]*math.cos(rotation) - v[0]*math.sin(rotation)))))
            
            colors.append((255,0,0))
            
def travel(verts, times, colors, start, end, duration = 0.035, steps = 4):
    start_time = times[-1]
    t = start_time
    
    for i in range(steps):
        t += duration / steps
        times.append(t)
        
        interp = (t - start_time) / duration
        verts.append(((1-interp)*start[0] + (interp)*end[0], (1-interp)*start[1] + (interp)*end[1]))
        
        colors.append((0,0,0))

def main():
    verts = [(0,0)]
    colors = [(0,0,0)]
    times = [0]
    
    # first green shot
    shootDisruptor(verts, times, colors, (250,15), (125,109))
    littleExplosion(verts, times, colors, (125,109))
    travel(verts, times, colors, (125,109), (250,15))
    
    # second green shot
    shootDisruptor(verts, times, colors, (250,15), (125,109))
    littleExplosion(verts, times, colors, (125,109))
    travel(verts, times, colors, (125,109), (78, 171))
    
    #phaser
    shootPhaser(verts, times, colors, (78, 171), (90,230), (140,200))
    travel(verts, times, colors, (140,200), (150, 100))
    travel(verts, times, colors, (150, 100), (250,15))
    
    
    # first green shot
    shootDisruptor(verts, times, colors, (250,15), (125,109))
    littleExplosion(verts, times, colors, (125,109))
    travel(verts, times, colors, (125,109), (250,15))
    
    # second green shot
    shootDisruptor(verts, times, colors, (250,15), (125,109))
    littleExplosion(verts, times, colors, (125,109))
    travel(verts, times, colors, (125,109), (250,15))
    
    # third green shot
    shootDisruptor(verts, times, colors, (250,15), (125,109))
    littleExplosion(verts, times, colors, (125,109))
    travel(verts, times, colors, (125,109), (0, 16), duration = 0.21)
    
    #shoot blue photon torpoedo
    #first point
    times.append(times[-1] + 0.02)
    verts.append((0, 16))
    colors.append((0,0,255))
    
    #second point
    times.append(times[-1] + 0.84)
    verts.append((75,75))
    colors.append((0,0,255))
    
    #third point
    times.append(times[-1] + 0.497)
    verts.append((150,25))
    colors.append((0,0,255))
    
    #fourth point
    times.append(times[-1] + 0.2)
    verts.append((188,109))
    colors.append((0,0,255))
    
    #big explosion
    bigExplosion(verts, times, colors, (188,78))
    travel(verts, times, colors, (225,96), (250, 15), duration = 0.0315)
    
    #last push to make sure it has rested
    times.append(times[-1] + (0.0315/4.0))
    verts.append((250,15))
    colors.append((0,0,0))
    
    interpolateAndOutput(verts[1:], colors[1:], times[1:])
    
if __name__ == "__main__":
    main()