from vpython import *

sf_sun = 10     
sf_earth = 300  
sf_mars = 250
t = 0
dt = 60*60
G = 6.673e-11
rate_speed = 300
dday = 10

R_sun = 696000000
R_earth = 6400000 
R_mars = 3389500

r_earth = 1.5e11
r_mars = 2.28e11

gap_earth = 0.8
gap_mars = 0.6

scene = canvas(title="Kepler's 2nd & 3rd Law Simulation", width=800, height=600)
scene.center = vec(0, 0, 0)

SUN = sphere(canvas=scene, pos=vector(0,0,0), radius=sf_sun*R_sun, color=color.orange)
SUN.mass = 1.99e30

EARTH = sphere(canvas=scene, pos=vector(r_earth,0,0), radius=sf_earth*R_earth, color=color.white, make_trail=True)
EARTH.mass = 5.97e24
vi_earth = sqrt(G*SUN.mass/r_earth)
EARTH.v = vec(0, vi_earth*gap_earth, 0)

MARS = sphere(canvas=scene, pos=vector(r_mars,0,0), radius=sf_mars*R_mars, color=color.red, make_trail=True, trail_color=color.red)
MARS.mass = 6.42e23
vi_mars = sqrt(G*SUN.mass/r_mars)
MARS.v = vec(0, vi_mars*gap_mars, 0)

SUN.v = -(EARTH.v*EARTH.mass + MARS.v*MARS.mass)/SUN.mass

scene.append_to_caption(f"\n[실시간 {dday}일 간격 면적 및 공전 데이터 로그]\n")
output_text = wtext(text="")
log_list = []
max_logs = 10

earth_polygons = []
mars_polygons = []

scene.waitfor('click')

old_pos_e = vec(EARTH.pos.x, EARTH.pos.y, EARTH.pos.z)
total_arc_area_e = 0

old_pos_m = vec(MARS.pos.x, MARS.pos.y, MARS.pos.z)
total_arc_area_m = 0

day_counter = 0

earth_max = r_earth
earth_min = r_earth
earth_old_y = EARTH.pos.y

mars_max = r_mars
mars_min = r_mars
mars_old_y = MARS.pos.y

earth_done = False
mars_done = False
k_e_saved = ""
k_m_saved = ""

while t < 10*365*24*60*60:
    rate(rate_speed)
    
    old_dist_e = mag(old_pos_e)
    rel_e = EARTH.pos - SUN.pos
    dist_e = mag(rel_e)
    
    if dist_e > earth_max: earth_max = dist_e
    if dist_e < earth_min: earth_min = dist_e
        
    EARTH.f = -G*SUN.mass*EARTH.mass/dist_e**2*norm(rel_e)
    EARTH.v = EARTH.v + EARTH.f/EARTH.mass*dt
    EARTH.pos = EARTH.pos + EARTH.v*dt

    dot_e = dot(old_pos_e, rel_e)
    mag_e = mag(old_pos_e) * mag(rel_e)
    cos_e = dot_e / mag_e
    if cos_e > 1.0: cos_e = 1.0
    elif cos_e < -1.0: cos_e = -1.0
    theta_e = acos(cos_e)
    
    avg_dist_e = (old_dist_e + dist_e) / 2
    arc_area_e = 0.5 * (avg_dist_e**2) * theta_e
    total_arc_area_e += arc_area_e

    ep0 = vec(0, 0, 0)
    ep1 = norm(old_pos_e) * avg_dist_e
    ep2 = norm(rel_e) * avg_dist_e
    ev0 = vertex(pos=ep0, color=color.yellow, opacity=0.3)
    ev1 = vertex(pos=ep1, color=color.yellow, opacity=0.3)
    ev2 = vertex(pos=ep2, color=color.yellow, opacity=0.3)
    tri_e = triangle(canvas=scene, v0=ev0, v1=ev1, v2=ev2)
    earth_polygons.append(tri_e)

    old_dist_m = mag(old_pos_m)
    rel_m = MARS.pos - SUN.pos
    dist_m = mag(rel_m)
    
    if dist_m > mars_max: mars_max = dist_m
    if dist_m < mars_min: mars_min = dist_m
        
    MARS.f = -G*SUN.mass*MARS.mass/dist_m**2*norm(rel_m)
    MARS.v = MARS.v + MARS.f/MARS.mass*dt
    MARS.pos = MARS.pos + MARS.v*dt

    dot_m = dot(old_pos_m, rel_m)
    mag_m = mag(old_pos_m) * mag(rel_m)
    cos_m = dot_m / mag_m
    if cos_m > 1.0: cos_m = 1.0
    elif cos_m < -1.0: cos_m = -1.0
    theta_m = acos(cos_m)
    
    avg_dist_m = (old_dist_m + dist_m) / 2
    arc_area_m = 0.5 * (avg_dist_m**2) * theta_m
    total_arc_area_m += arc_area_m

    mp0 = vec(0, 0, 0)
    mp1 = norm(old_pos_m) * avg_dist_m
    mp2 = norm(rel_m) * avg_dist_m
    mv0 = vertex(pos=mp0, color=color.cyan, opacity=0.3)
    mv1 = vertex(pos=mp1, color=color.cyan, opacity=0.3)
    mv2 = vertex(pos=mp2, color=color.cyan, opacity=0.3)
    tri_m = triangle(canvas=scene, v0=mv0, v1=mv1, v2=mv2)
    mars_polygons.append(tri_m)

    SUN.pos = SUN.pos + SUN.v*dt

    if not earth_done and earth_old_y < 0 and EARTH.pos.y >= 0:
        a_e = (earth_max + earth_min) / 2
        k_e = (t**2) / (a_e**3)
        k_e_saved = f" [지구 공전완료! T^2/a^3: {k_e:.4e}, T={t}, a={a_e}]"
        earth_done = True
        
    if not mars_done and mars_old_y < 0 and MARS.pos.y >= 0:
        a_m = (mars_max + mars_min) / 2
        k_m = (t**2) / (a_m**3)
        k_m_saved = f" [화성 공전완료! T^2/a^3: {k_m:.4e}, T={t}, a={a_m}]"
        mars_done = True

    old_pos_e = vec(rel_e.x, rel_e.y, rel_e.z)
    old_pos_m = vec(rel_m.x, rel_m.y, rel_m.z)

    if t % (24*60*60) == 0 and t > 0:
        day_counter += 1
        
        if day_counter == dday:
            current_day = int(t / (24*60*60))
            past_day = current_day - dday
            log_line = f"{past_day}~{current_day}일째 면적 ({dday}일 누적) -> 지구: {total_arc_area_e:.8e} | 화성: {total_arc_area_m:.8e} |"
            
            if k_e_saved:
                log_line += k_e_saved
            if k_m_saved:
                log_line += k_m_saved
                
            log_list.insert(0, log_line)
            
            if len(log_list) > max_logs:
                log_list.pop()
                
            output_text.text = "\n".join(log_list)
            total_arc_area_e = 0
            total_arc_area_m = 0
            day_counter = 0
            
            for p in earth_polygons: p.visible = False
            for p in mars_polygons: p.visible = False
            earth_polygons.clear()
            mars_polygons.clear()

    earth_old_y = EARTH.pos.y
    mars_old_y = MARS.pos.y
    
    t = t + dt
