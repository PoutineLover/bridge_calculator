import copy
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

BRIDGE_LENGTH = 1290
class EntireBridge: 

    # Instance Variables
    # cross_sections: List of cross-sections that make up the bridge
    # length: Length of the bridge
    # d_locations: List of indices where a diaphragm exists
    # initial_members: Basic top, bottom and side members
    def __init__(self, cross_section_members, d_spacing, length):
        # Number of cross sections of the bridge
        self.cross_section = CrossSection(cross_section_members)
        # Indices where a diaphragm exists
        self.d_spacing = d_spacing
        # Length of the bridge in mm
        self.length = length

class CrossSection:

    # Instance Variables
    # members_dict: dictionary of Beam objects representing the cross-section
    # y_bar: y_bar
    # I: I
    def __init__(self, members):
        self.members_dict = copy.deepcopy(members)

        # If there is a diaphragm at this spot, it adds it as a separate member between
        # the left and right members of the cross-section

        self.height = self.get_height()
        self.y_bar = self.get_y()
        self.I = self.get_I()

    # Returns the height of the cross_section
    def get_height(self):
        return self.members_dict["t"].y_pos + self.members_dict["t"].h/2

    # Output: Returns the y_bar of the "members_dict"
    def get_y(self): 
        total_ad, total_area = 0, 0
        for beam in self.members_dict.values():
            total_ad += beam.y_pos * beam.b * beam.h
            total_area += beam.b * beam.h
            #print(f'New ad: {beam.y_pos * beam.b * beam.h}')
            #print(f'New area: {beam.b * beam.h}')
        return total_ad/total_area

    # Output: returns I of the "members_dict"
    def get_I(self): # Return list
        y_bar = self.y_bar
        total_I = 0
        for beam in self.members_dict.values():
            total_I += (beam.b * beam.h**3)/12 + (beam.b * beam.h) * (beam.y_pos - y_bar)**2 # I_0 + Ad^2
        return total_I

    # Returns the Q value at a specific height
    def get_q_at_height(self, height):
        member_list_copy = copy.deepcopy(self.members_dict) # Deep copy needed bc list holds other objects

        # The list of members is adjusted to only include portions of beams 
        # above height
        relevant_beams = dict(filter(lambda el: el[1].y_pos + el[1].h/2 > height, member_list_copy.items()))
        to_list = list(relevant_beams.values())
        
        total_Q = 0
        # Calculates Q, accounting for parts of each beam that extend below "height"
        for beam in to_list:
            top_y = (beam.y_pos + beam.h/2) # y coordinate of topmost point of each beam beam
            bottom_y = max(beam.y_pos - beam.h/2, height) # y coordinate of bottommost point of new beam

            beam.y_pos = (bottom_y + top_y)/2 # New y coordinate of beam
            beam.h = top_y - bottom_y # New height of beam
            # Adds A * d
            total_Q += beam.b * beam.h * (beam.y_pos - self.y_bar)
        return total_Q

    # Gets the total width at a height h of the cross section
    # Useful for shear calculations
    # If type = "+", the height is shifted 1E-06 up
    # If type = "-", the height is shifted 1E-06 down
    def get_b_at_height(self, height, type):
        if type == '+':
            height = height + 0.001
        elif type == '-':
            height = height - 0.001
        # Finds beams in the cross section whose area intersects y = height
        relevant_beams = []
        for beam in self.members_dict.values():
            if beam.y_pos + beam.h/2 >= height and beam.y_pos - beam.h/2 <= height:
                relevant_beams.append(beam)    
        #relevant_beams = dict(filter(lambda el: el[1].y_pos + el[1].h/2 >= height and el[1].y_pos - el[1].h/2 <= height, self.members_dict))
        total_b = 0
        for beam in relevant_beams:
            total_b += beam.b
        return total_b

    def to_string(self):
        scientific_notation = "{:.4e}".format(self.I)
        scientific_notation_2 = "{:.4e}".format(self.get_q_at_height(self.y_bar))
        scientific_notation_3 = "{:.4e}".format(self.get_q_at_height(self.members_dict['t'].y_pos - (self.members_dict['t'].h/2)))
        b1 = self.get_b_at_height(self.y_bar, "0")
        print(f'y_bar: {self.y_bar} mm')
        print(f'I: {scientific_notation} mm^4')
        print(f'Q @ y_bar: {scientific_notation_2}')
        print(f'Q @ top glue area: {scientific_notation_3}')
        print(f"Base at y_bar: {b1}")
        #print(f"Height at glue fold: {b2}")

    def draw_cross_section(self):
        plt.axes()
        for beam in self.members_dict.values():
            rectangle = plt.Rectangle((beam.x_pos - beam.b/2, beam.y_pos - beam.h/2), beam.b, beam.h, fc="none", ec="red", linewidth=1)
            plt.gca().add_patch(rectangle)
            plt.axis('scaled')
        plt.show()

class Beam:
    # x is the x-coordinate of the center of the beam
    # y is the y-coordinate of ""
    # b is the width of the beam
    # h is the height of the beam
    def __init__(self, x, y, b, h):
        self.x_pos = x
        self.y_pos = y
        self.b = b
        self.h = h

    # Returns a copy of the beam
    def copy(self):
        return Beam(self.x_pos, self.y_pos, self.b, self.h)

    def to_string(self):
        print(f"x_pos: {self.x_pos}\ny_pos: {self.y_pos}\nwidth: {self.b}\nheight: {self.h}\n")

# ----- SHEAR BUCKLING ----- #

# Gets the smallest shear buckling force needed to break the bridge
def VFailBuck(bridge, E, mu):
    cross_section = copy.deepcopy(bridge.cross_section)
    diaphragm_list = copy.deepcopy(bridge.d_spacing)
    # For loop that finds the maximum distance between two diaphragms
    max_distance = -math.inf
    for i in range(1, len(diaphragm_list)):
        distance = diaphragm_list[i] - diaphragm_list[i - 1]
        max_distance = max(max_distance, distance)
    # Height of the bridge
    height = cross_section.height 
    # Calculate
    shear_crit = ((5 * math.pi**2 * E) / (12 * (1 - mu**2))) * ((1.27/max_distance)**2 + (1.27/height)**2)

    I = cross_section.I
    y_bar = cross_section.y_bar
    b = cross_section.get_b_at_height(y_bar, "0")
    Q = cross_section.get_q_at_height(y_bar)
    V = shear_crit * I * b / Q
    return V

    # ----- OTHER SHEAR FAILURES -----#
    # Returns a list of 3 elements
# They correspond to shear failure of material, shear failure of top glue fold, and
# shear failure of bottom glue fold
def V_fail(bridge, tau_u, tau_g):

    # V_fail_list is a list of length 3 that contains V_fail_mat and V_fail_glue (top and bot)
    V_fail_list = []
    cross_section = copy.deepcopy(bridge.cross_section)

    members_dict = cross_section.members_dict
    I = cross_section.I
    y_bar = cross_section.y_bar
    total_height = cross_section.height

    #------- Shear Material Failure -------#
    Q = cross_section.get_q_at_height(y_bar)
    b = cross_section.get_b_at_height(y_bar, '0')
    V = tau_u * I * b / Q
    V_fail_list.append(V)

    #------- Shear Glue Failure At TOP -------#

    #------- Connection from fold to top flange -------#
    height = total_height - members_dict['t'].h
    Q = cross_section.get_q_at_height(height)
    b = cross_section.get_b_at_height(height, '-') - (1.27 * 2)
    V = tau_g * I * b / Q
    V_fail_list.append(V)

    #------- Shear Glue Failure At BOTTOM -------#

    #------- Connection from fold to bottom flange -------#
    # Checks to see if a bottom fold even exists
    if not members_dict.get("f3", 0) == 0:
        height = members_dict['b'].h
        Q = cross_section.get_q_at_height(height)
        b = cross_section.get_b_at_height(height, '+') - (1.27 * 2)
        V = tau_g * I * b / Q
        V_fail_list.append(V)

    return V_fail_list

    # ----- PLATE BUCKLING ----- #
# Returns a list containing 3 lists, each sub-list represents Cases 1-3 for plate buckling
# The 1st element in each sublist is the minimum M for when M < 0
# The 2nd element in each sublist is the minimum M for when M > 0
# Note that both elements are still positive
def M_Fail_Buck(bridge, E, mu, BMD):
    # List of length 3 where each row holds n values of M_fail for either case 1, case 2, or case 3 
    # Row 0 -> Case 1
    # Row 1 -> Case 2
    # Row 2 -> Case 3
    M_buck_list = []
    # Case 1, buckling with fixed ends
    M_buck_c1 = []
    # Case 2, buckling with one free end
    M_buck_c2 = []
    # Case 3, buckling with fixed ends, linearly varying comp. stress
    M_buck_c3 = []
    
    cross_section = copy.deepcopy(bridge.cross_section)
    # Used in all case calculations
    members_dict = cross_section.members_dict
    height = cross_section.height
    y_bar = cross_section.y_bar
    I = cross_section.I

    #----------- CASE 1 ------------#

    # If the BMD < 0, then the bottom flange is in compression
    t1 = members_dict['b'].h
    b1 = abs(members_dict['r'].x_pos - members_dict['l'].x_pos)
    y =  y_bar
    sigma_crit1 = ((4 * (math.pi ** 2) * E) / (12 * (1 - mu ** 2))) * ((t1 / b1) ** 2)
    M = (sigma_crit1 * I) / y
    M_buck_c1.append(M)

    # If the BMD > 0, then the top flange is in compression
    t1 = members_dict['t'].h
    b1 = abs(members_dict['r'].x_pos - members_dict['l'].x_pos)
    y = height - y_bar
    sigma_crit1 = ((4 * (math.pi ** 2) * E) / (12 * (1 - mu ** 2))) * ((t1 / b1) ** 2)
    M = (sigma_crit1 * I) / y
    M_buck_c1.append(M)
    
    #----------- CASE 2 ------------#

    # If the BMD < 0 , then the bottom flange is in comopression
    t2 = members_dict['b'].h
    b2 = (members_dict['b'].b - (members_dict['r'].x_pos - members_dict['l'].x_pos) / 2)
    y = y_bar
    sigma_crit2 = ((0.425 * (math.pi ** 2) * E) / (12 * (1 - mu ** 2))) * ((t2 / b2) ** 2)
    M = (sigma_crit2 * I) / y
    M_buck_c2.append(M)

    # If the BMD > 0, then the top flange is in compression
    t2 = members_dict['t'].h
    b2 = (members_dict['t'].b - (members_dict['r'].x_pos - members_dict['l'].x_pos)) / 2 
    y = height - y_bar
    sigma_crit2 = ((0.425 * (math.pi ** 2) * E) / (12 * (1 - mu ** 2))) * ((t2 / b2) ** 2)
    M = (sigma_crit2 * I) / y
    M_buck_c2.append(M)
            
    #----------- CASE 3 ------------#

    # If BMD < 0, then the bottom portion of the web is in compression
    t3 = members_dict['l'].b
    b3 = y_bar - (members_dict['l'].y_pos - members_dict['l'].h/2)
    y = b3
    sigma_crit3 = ((6 * (math.pi ** 2) * E) / (12 * (1 - mu**2))) * ((t3 / b3) ** 2)
    M = sigma_crit3 * I / y
    M_buck_c3.append(M) 

    # If BMD > 0, then the top portion of the web is in compression
    t3 = members_dict['l'].b
    b3 = (members_dict['t'].y_pos - (members_dict['t'].h)/2) - y_bar
    y = b3
    sigma_crit3 = ((6 * (math.pi ** 2) * E) / (12 * (1 - mu**2))) * ((t3 / b3) ** 2)
    M = sigma_crit3 * I / y
    M_buck_c3.append(M) 
                
    M_buck_list.append(M_buck_c1)
    M_buck_list.append(M_buck_c2)
    M_buck_list.append(M_buck_c3)
    return M_buck_list

# ----- MOMENT FLEXURAL FAILURES ----- #

# Returns the minimum M needed to break the bridge when M < 0 and when M > 0
# These produce different values because the top and the bottom break at different moments
# Output: A tuple of length 2; 
# The first element is the minimum M needed to break the bridge when M < 0
# The second element is the minimum M needed to break the bridge when M > 0
def MfailMatT(bridge, BMD, sigT):
    M_Fail_T_list = []
    cross_section = copy.deepcopy(bridge.cross_section)
    I = cross_section.I
    y_bar = cross_section.y_bar
    height = cross_section.height
    # When BMD < 0, the top of the beam is in tensile stress
    # Tensile Stress for the top of the beam
    M_Fail_T_list.append(sigT * I / (height - y_bar))
    # When BMD > 0, the bottom of the beam is in tensile stress
    # Tensile Stress for the bottom of the beam
    M_Fail_T_list.append(sigT * I / y_bar)
    return M_Fail_T_list

def MfailMatC(bridge, BMD, sigC):
    M_Fail_C_list = []
    cross_section = copy.deepcopy(bridge.cross_section)
    I = cross_section.I
    y_bar = cross_section.y_bar
    height = cross_section.height

    # When BMD < 0, the bottom of the beam is in compressive stress
    M_Fail_C_list.append(sigC * I / y_bar)
    # When BMD > 0, the top of the beam is in compressive stress
    M_Fail_C_list.append(sigC * I / (height - y_bar))
    return M_Fail_C_list

# ----- CREATING THE TRAIN PLACEMENTS ----- #
def get_train_placements(train_loads_list, length, step):
    train_list = copy.deepcopy(train_loads_list)
    # Shifts the train back so that the rightmost wheel is at position 0
    for i in range(0, len(train_list)):
        train_list[i] -= train_list[len(train_list) - 1]

    # For loop that runs from 0 to the length of the bridge + position of right most wheel
    # Kind of HARD CODED
    LENGTH = length + train_loads_list[-1]
    STEP = step
    NUM_ITERATIONS = math.ceil(LENGTH/STEP)
    list_of_train_placements = []
    for i in range(0, NUM_ITERATIONS):
        for index in range(0, len(train_list)):
            train_list[index] += STEP
        # Only keeps loads that are on the train (HARD CODED)
        train_list_loop = list(filter(lambda el: el >= 0 and el <= BRIDGE_LENGTH, train_list))
        # Makes the SF_points list and appends it to the list of train_placements
        sf_points = []
        for position in train_list_loop:
            sf_points = make_SF_points(position, 400 / 6, sf_points)
        list_of_train_placements.append(sf_points)
        #sf_tuples = SFDBMD.get_SF_tuples(sf_points)
        #bm_tuples = SFDBMD.get_BM_tuples(sf_tuples)
        #SFDBMD.draw_SFD(sf_tuples)
        #SFDBMD.draw_BMD(bm_tuples)
    return list_of_train_placements

# ----- CREATING THE SFD and BMDS ----- #
# Output: A list of SF_points which convey the following:
# Positions on the bridge where a point force exists. Each element of SF_points is a tuple,
# whose first element is a location on the bridge, and whose second element is the point force 
# being applied
# Positive values are Ay and By
# Negative values are point loads being applied to the bridge
# The bridge will be positioned similar to how Design0 is positioned
def make_SF_points(xP, P, SF):
    
    # Taking the moment about point A
    # Adjust this so that it is no longer hard coded?
    By = P * (xP - 15) / 1060
    Ay = P - By

    #print(Ay)
    #print(By)

    # This list contains tuples of length 2 whose first element represents the x location 
    # of the shear force, and whose second element represents its shear force
    SF_points = copy.deepcopy(SF)

    # If the list is empty, the values are simply added in
    if len(SF_points) == 0:
        SF_points.append([15, Ay])
        SF_points.append([1075, By])
        SF_points.append([xP, -P])
    # Otherwise, the values are added to the SFD
    else:
        # Ay and By are always the first two forces on the list
        SF_points[0][1] += Ay
        SF_points[1][1] += By
        matched = False
        # Checks to see if there already exists a point load at xP
        for i in range(2, len(SF_points)):
            if xP == SF_points[i][0]:
                SF_points[i][1] += -P
                matched = True
                break
        # If not, it is given its own spot on the list
        if not matched:
            SF_points.append([xP, -P])

        ### ------ MAKING THE BMD FROM SF_points ----- ###

    return SF_points

# Takes SF_points generated by SFDBMD
# Returns a list of tuples of length 2. About each tuple:
# The first value is a location on the bridge
# The second value is the shear force present from that location to the location of the tuple
# before it (or zero if it is the first tuple)
def get_SF_tuples(SF_points):
    SF_list = copy.deepcopy(SF_points)
    sorted_SF_list = sorted(SF_list, key=lambda f: f[0], reverse=False)
    sorted_SF_list.append([BRIDGE_LENGTH, 0])
    sorted_SF_list.insert(0, [0, 0])
    # In a new list, each value is now converted to the sum of the forces that appear before it
    SF_tuples = []
    for i in range(0, len(sorted_SF_list)):
        total_shear = 0
        for j in range(0, i):
            total_shear += sorted_SF_list[j][1]
        SF_tuples.append((sorted_SF_list[i][0], total_shear))
    return SF_tuples

# Input: The list of tuples generated by get_SF_tuples
# Output: A list of tuples of length 2. About each tuple:
# The first value is a location along the length of the bridge
# The second value is the Bending Moment at that position
# Note that location of each point in the tuple is at a "notable" Bending Moment i.e there is a corner
def get_BM_tuples(SF_list_of_tuples):
    SF_tuples = copy.deepcopy(SF_list_of_tuples)
    BM_tuples = []
    total_moment = 0
    for i in range(0, len(SF_tuples)):
        # A (0, 0) tuple is made at the beginning
        if i == 0:
            BM_tuples.append((0, 0))
        # Calculates the area of the shear sections and 
        else:
            # total_moment += base * height of shear rectangle section
            total_moment += (SF_tuples[i][0] - SF_tuples[i - 1][0]) * SF_tuples[i][1]
            # Gets rid of annoying floating points
            if abs(total_moment) < 0.001:
                total_moment = 0
            BM_tuples.append((SF_tuples[i][0], total_moment))
    #print(BM_tuples)
    return BM_tuples

# Returns the points in the BMD where a zero occurs. This is important if 
# we need to plot the BMD for various failures. 
# Output: A list of tuples of length 2
# The first element of each tuple is the x-coordinate of the zero
# The second element is the sign change:
# 1 if the trend is neg -> pos
# 0 if the trend is pos -> neg
def find_bm_zeroes(BM_list_of_tuples):
    list_of_zeroes = []
    for i in range(0, len(BM_list_of_tuples) - 1):
        # Checking if a sign change took place
        if (BM_list_of_tuples[i][1] >= 0 and BM_list_of_tuples[i + 1][1] < 0) or (BM_list_of_tuples[i][1] < 0 and BM_list_of_tuples[i + 1][1] < 0):
        # Does linear interpolation to find the zero
            m = (BM_list_of_tuples[i + 1][1] - BM_list_of_tuples[i][1]) / (BM_list_of_tuples[i + 1][0] - BM_list_of_tuples[i][0])
            # Process of finding the zero:
            # y_2 - y_1 = m(x_2 - x_1)
            # 0 - y_1 = m(x - x_1)
            # mx = mx_1 - y_1
            # x = (mx_1 - y_1)/m
            x = (m * BM_list_of_tuples[i][0] - BM_list_of_tuples[i][1]) / m
            # neg to pos
            if m > 0:
                list_of_zeroes.append((x, 1))
            else:
                list_of_zeroes.append((x, 0))
    return list_of_zeroes

# Takes in the list of tuples generated by get_SF_tuples and draws an SFD using those points
def draw_SFD(SF_list_of_tuples):
    # The goal is to have every "corner" of the SFD as a point in the list.
    # The following lines of code add these missing corners
    SF_tuples = copy.deepcopy(SF_list_of_tuples)
    SF_tuples.reverse()
    #SF_tuples.append([0, 0])
    #print(SF_tuples)

    # Formats the sorted SFD list into a list of important "corners" 
    for i in range(len(SF_tuples) - 2, -1, -1):
        SF_tuples.insert(i + 1, (SF_tuples[i + 1][0], SF_tuples[i][1]))
    #print(SF_tuples)
    
    # Draws a zero axis in black
    plt.plot([0, BRIDGE_LENGTH], [0, 0], 'k--')
    # Draws the SFD using plt
    for i in range(len(SF_tuples) - 1, 0, -1):
        x_coordinates = [SF_tuples[i][0], SF_tuples[i - 1][0]]
        y_coordinates = [SF_tuples[i][1], SF_tuples[i - 1][1]]
        plt.plot(x_coordinates, y_coordinates, 'r')
    # Adds a legend for the SFD
    sfd_patch = patches.Patch(color='red', label='SFD Of Bridge')
    sfd_line = [Line2D([0], [0], color='red', linewidth=1, linestyle='-')]
    sfd_label = ["SFD Of Bridge"]
    #plt.legend(sfd_line, sfd_label)
    #plt.show()

# Takes in the list of tuples generated by get_BM_tuples and draws a BMD using those points
def draw_BMD(BM_list_of_tuples):
    BM_tuples = copy.deepcopy(BM_list_of_tuples)
    # Unlike SFD_draw, we don't need to add additional points
    plt.plot([0, BRIDGE_LENGTH], [0, 0], 'k--')
    for i in range(0, len(BM_tuples) - 1):
        x_coordinates = [BM_tuples[i][0], BM_tuples[i + 1][0]]
        y_coordinates = [BM_tuples[i][1], BM_tuples[i + 1][1]]
        plt.plot(x_coordinates, y_coordinates, 'b')
    bmd_line = [Line2D([0], [0], color='blue', linewidth=1, linestyle='-')]
    bmd_label = ["BMD Of Bridge"]
    #plt.legend(bmd_line, bmd_label)

# Checks whether the BMD at this point is positive or negative
# It does this by finding the two positions of BMD that the position is 
# in between, and then does a linear approximation (y = mx + b)
# Note that BM_list_of_tuples is a list of tuples containing the 
# Bending Moment at maximum points and endpoints, not a list of values
def check_BMD_signage(position, BM_list_of_tuples):
    BM_tuples = copy.deepcopy(BM_list_of_tuples)
    #print(BM_tuples)
    for i in range(0, len(BM_tuples) - 1):
        # When the two positions in the BMD have been found
        if position >= BM_tuples[i][0] and position <= BM_tuples[i + 1][0]:
            x_coordinates = (BM_tuples[i][0], BM_tuples[i + 1][0])
            y_coordinates = (BM_tuples[i][1], BM_tuples[i + 1][1])
            #print(x_coordinates)
            #print(y_coordinates)
            # If delta-x is zero, set the distance equal to 0.001
            if (x_coordinates[1] - x_coordinates[0]) < 0.001:
                m = (y_coordinates[1] - y_coordinates[0]) / 0.001
            else:
                m = (y_coordinates[1] - y_coordinates[0]) / (x_coordinates[1] - x_coordinates[0])
            # Gets the distance between position and left endpoint and moves that amount * slope
            # from left endpoint
            distance = position - x_coordinates[0]
            if y_coordinates[0] + m * distance > 0:
                return 1
                break
            else: 
                return -1
                break
    return "Nothing!"

# Gets the maximum absolute shear force present in sf_tuples
# Note: Returns the full tuple
def get_max_SF(sf_tuples):
    return max(sf_tuples, key = lambda s_force: abs(s_force[1]))

# Gets the maximum or minimum bending moment in bm_tuples
# Note: Returns the full tuple
def get_max_BM(bm_tuples, type):
    if type == "p":
        return max(bm_tuples, key = lambda b_moment: b_moment[1])
    elif type == "n":
        return min(bm_tuples, key = lambda b_moment: b_moment[1])
    else:
        return "BAD INPUT"

def find_lowest_fos(bridge, train_loads_list, bridge_length):
    # Constants
    SigT = 30
    SigC = 6
    E = 4000
    TauU = 4
    TauG = 2
    mu = 0.2

    # Getting shear failures
    fail_shear = V_fail(bridge, TauU, TauG)
    fail_v_mat = fail_shear[0]
    fail_v_glue_top = fail_shear[1]
    fail_v_buck = VFailBuck(bridge, E, mu)
    min_v_fail = min(fail_v_mat, fail_v_glue_top, fail_v_buck)
    if len(fail_shear) > 2:
        fail_v_glue_bot = fail_shear[2]
        min_v_fail = min(fail_v_mat, fail_v_glue_top, fail_v_glue_bot, fail_v_buck)

    # Getting moment failures
    fail_m_comp = MfailMatC(bridge, None, SigC)
    fail_m_tens = MfailMatC(bridge, None, SigT)
    fail_m_buck_1 = M_Fail_Buck(bridge, E, mu, None)[0]
    fail_m_buck_2 = M_Fail_Buck(bridge, E, mu, None)[1]
    fail_m_buck_3 = M_Fail_Buck(bridge, E, mu, None)[2]
    min_m_fail_neg = min(fail_m_comp[0], fail_m_tens[0], fail_m_buck_1[0], fail_m_buck_2[0], fail_m_buck_3[0])
    min_m_fail_pos = min(fail_m_comp[1], fail_m_tens[1], fail_m_buck_1[1], fail_m_buck_2[1], fail_m_buck_3[1])

    list_of_train_loads = get_train_placements(train_loads_list, bridge_length, 10)
    min_v_fos = math.inf
    min_m_pos_fos = math.inf
    min_m_neg_fos = math.inf
    for sf_configuration in list_of_train_loads:
        # Makes the SFD and BMD configurations
        sf_tuples = get_SF_tuples(sf_configuration)
        bm_tuples = get_BM_tuples(sf_tuples)
        # Gets the absolute maximum shear force and largest positive and negative
        # bending moments present in the bridge
        max_sf = get_max_SF(sf_tuples)[1]
        max_pos_bm = get_max_BM(bm_tuples, "p")[1]
        max_neg_bm = get_max_BM(bm_tuples, "n")[1]

        # If the maximum force is 0, it won't try dividing by 0
        v_fos = math.inf
        m_fos_pos = math.inf
        m_fos_neg = math.inf
        if not max_sf == 0:
            v_fos = min_v_fail/abs(max_sf)
        if not max_pos_bm == 0:
            m_fos_pos = min_m_fail_pos/abs(max_pos_bm)
        if not max_neg_bm == 0:
            m_fos_neg = min_m_fail_neg/abs(max_neg_bm)

        min_v_fos = min(min_v_fos, v_fos)
        min_m_pos_fos = min(min_m_pos_fos, m_fos_pos)
        min_m_neg_fos = min(min_m_neg_fos, m_fos_neg)
        
    print(f'Minimum V FOS: {min_v_fos}')
    print(f'Minimum positive M FOS: {min_m_pos_fos}')
    print(f'Minimum negative M FOS: {min_m_neg_fos}')

# ----- FINDING THE MINIMUM POINT LOAD THAT WILL BREAK THE BRIDGE ----- #
def find_max_P(bridge):
    # Constants
    SigT = 30
    SigC = 6
    E = 4000
    TauU = 4
    TauG = 2
    mu = 0.2

    # Getting shear failures
    fail_v_buck = VFailBuck(bridge, E, mu)
    fail_shear = V_fail(bridge, TauU, TauG)
    fail_v_mat = fail_shear[0]
    fail_v_glue_top = fail_shear[1]
    min_v_fail = min(fail_v_mat, fail_v_glue_top, fail_v_buck)
    if len(fail_shear) > 2:
        fail_v_glue_bot = fail_shear[2]
        min_v_fail = min(min_v_fail, fail_v_glue_bot)

    # Getting moment failures
    fail_m_comp = MfailMatC(bridge, None, SigC)
    fail_m_tens = MfailMatC(bridge, None, SigT)
    fail_m_buck_1 = M_Fail_Buck(bridge, E, mu, None)[0]
    fail_m_buck_2 = M_Fail_Buck(bridge, E, mu, None)[1]
    fail_m_buck_3 = M_Fail_Buck(bridge, E, mu, None)[2]
    min_m_fail_neg = min(fail_m_comp[0], fail_m_tens[0], fail_m_buck_1[0], fail_m_buck_2[0], fail_m_buck_3[0])
    min_m_fail_pos = min(fail_m_comp[1], fail_m_tens[1], fail_m_buck_1[1], fail_m_buck_2[1], fail_m_buck_3[1])
    print(min_m_fail_neg)
    print(min_m_fail_pos)

    # Initial P Load of 
    sf_points = make_SF_points(565, 1, [])
    sf_points = make_SF_points(bridge.length - 15, 1, sf_points)
    sf_tuples = []
    bm_tuples = []
    counter = 2
    while True:
        sf_tuples = get_SF_tuples(sf_points)
        bm_tuples = get_BM_tuples(sf_tuples)

        # Getting maximum values of shear and bending moment
        sf_max = get_max_SF(sf_tuples)
        bm_max_pos = get_max_BM(bm_tuples, "p") 
        bm_max_neg = get_max_BM(bm_tuples, "n") 
        #print(f'Max shear: {sf_max}')
        #print(f'Max pos M: {bm_max_pos}')
        #print(f'Max neg M: {bm_max_neg}')
    
        if sf_max[1] > min_v_fail:
            print("Shear Failure")
            break
        elif bm_max_pos[1] > min_m_fail_pos:
            print("Positive Moment Failure")
            break
        elif abs(bm_max_neg[1]) > min_m_fail_neg:
            print("Negative Moment Failure")
            break

        # Adds 1 more N to each point on the bridge
        sf_points = make_SF_points(565, 1, sf_points)
        sf_points = make_SF_points(bridge.length - 15, 1, sf_points)
        counter += 2
        
    print(f"Failed at: {counter}N")

    # ----- Drawing the SFD and BMD diagrams ----- #
def draw_other_BMDs(bm_tuples, main_bm, failure_tuples, bridge_length, labels, line_colors):
    bm_tuple_zeroes = copy.deepcopy(main_bm)
    # Scuffed fix for connecting the ends of the BMD. It assumes the BMD
    # ends going from positive to negative
    bm_tuple_zeroes.insert(0, (0, 0))
    bm_tuple_zeroes.append((bridge_length, 1))
    plt.plot([0, bridge_length], [0, 0], 'k--')
    # This nested for loop uses the bm_tuple_zeroes to help plot the moment failure diagrams
    for i in range(1, len(bm_tuple_zeroes)):
        # All of the moment piecewise functions are being drawn concurrently
        for index in range(0, len(failure_tuples)):
            x_coordinates = [bm_tuple_zeroes[i - 1][0], bm_tuple_zeroes[i][0]]
            # If the BMD is going upward, then that means it was negative before
            # and we want the negative moment value before that point
            # The signage of negative moments also needs to be physically set here
            if bm_tuple_zeroes[i][1] == 1:
                y_coordinates = [-failure_tuples[index][0], -failure_tuples[index][0]]
            # If the BMD is going downward, then that means it was positive before
            # and we want the positive moment value before that point
            else:
                y_coordinates = [failure_tuples[index][1], failure_tuples[index][1]]
            # Draws the horizontal line
            plt.plot(x_coordinates, y_coordinates, line_colors[index])
            # Draws the vertical line connecting the two horizontal lines
            x_coordinates = [bm_tuple_zeroes[i - 1][0], bm_tuple_zeroes[i - 1][0]]
            y_coordinates = [-failure_tuples[index][0], failure_tuples[index][1]]
            plt.plot(x_coordinates, y_coordinates, line_colors[index])
    # Draws the BMD
    draw_BMD(bm_tuples)
    # Makes the line symbols for the legend
    line_symbols = []
    for line_color in line_colors:
        line_symbols.append(Line2D([0], [0], color= line_color, linewidth=1, linestyle='-'))
    line_symbols.append(Line2D([0], [0], color= "b", linewidth=1, linestyle='-'))
    # Makes the labels for the legend
    legend_labels = labels
    legend_labels.append("BMD")
    plt.legend(line_symbols, legend_labels)

def draw_other_SFDs(sf_tuples, shear_failures, bridge_length, labels, line_colors):
    for i in range(0, len(shear_failures)):
        y_coordinates_pos = [shear_failures[i], shear_failures[i]]
        y_coordinates_neg = [-shear_failures[i], -shear_failures[i]]
        x_coordinates = [0, bridge_length]
        plt.plot(x_coordinates, y_coordinates_pos, line_colors[i])
        plt.plot(x_coordinates, y_coordinates_neg, line_colors[i])
    # Draws the SFD
    draw_SFD(sf_tuples)
    # Makes the line symbols for the legend
    line_symbols = []
    for line_color in line_colors:
        line_symbols.append(Line2D([0], [0], color= line_color, linewidth=1, linestyle='-'))
    line_symbols.append(Line2D([0], [0], color= "r", linewidth=1, linestyle='-'))
    # Makes the labels for the legend
    legend_labels = labels
    legend_labels.append("SFD")
    plt.legend(line_symbols, legend_labels)

if __name__ == "__main__":
    # Constants
    SigT = 30
    SigC = 6
    E = 4000
    TauU = 4
    TauG = 2
    mu = 0.2

    # LENGTH OF BRIDGE
    BRIDGE_LENGTH = 1290
    #------------------ DIAPHRAGM LOCATIONS (indices) --------------------# 
    d_spacing = [0, 184, 368, 553, 737, 921, 1289]

    #------------------ MEMBER DIMENSIONS --------------------#

    members_dict = {
        "t": Beam(0, 68.81, 100, 2.54), # Top
        "b": Beam(0, 1.27, 80, 2.54), # Bottom
        "l": Beam(-39.365, 35.04, 1.27, 65), # Left
        "r": Beam(39.365, 35.04, 1.27, 65), # Right
        "f1": Beam(33.73, 66.905, 10, 1.27), # Top Left Fold
        "f2": Beam(-33.73, 66.905, 10, 1.27), # Top Right Fold
        "f3": Beam(-33.73, 3.175, 10, 1.27), # Bottom Left Fold
        "f4": Beam(33.73, 3.175, 10, 1.27) # Bottom Left Fold
    }
    bridge = EntireBridge(members_dict, d_spacing, BRIDGE_LENGTH)

    # Train load list
    train_loads_list = [0, 176, 340, 516, 680, 856]

    # Get FOS when train is going across
    find_lowest_fos(bridge, train_loads_list, BRIDGE_LENGTH)

    # Get largest P that can go on the bridge (P is load of both points combined)
    find_max_P(bridge)

    # Getting the SFDS and BMDS of the bridge

    sf_points = make_SF_points(545, 475, [])
    sf_points = make_SF_points(BRIDGE_LENGTH - 15, 475, sf_points)
    sf_tuples = get_SF_tuples(sf_points)
    bm_tuples = get_BM_tuples(sf_tuples)
    bm_zeroes = find_bm_zeroes(bm_tuples)

    # Drawing the flexural stress failures
    mat_comp_fail = MfailMatC(bridge, None, SigC)
    mat_tens_fail = MfailMatT(bridge, None, SigT)
    flexural_failures = [mat_comp_fail, mat_tens_fail]
    flexural_failure_colors = ["chocolate", "green"]
    flexural_failure_labels = ["Compressive Flexural Stress", "Tensile Flexural Stress"]
    draw_other_BMDs(bm_tuples, bm_zeroes, flexural_failures, BRIDGE_LENGTH, flexural_failure_labels, flexural_failure_colors)
    plt.show()
    plt.clf()
    # Drawing the plate buckling stress failures
    moment_buckling = M_Fail_Buck(bridge, E, mu, None)
    buck_case_1 = moment_buckling[0]
    buck_case_2 = moment_buckling[1]
    buck_case_3 = moment_buckling[2]
    buckling_failures = [buck_case_1, buck_case_2, buck_case_3]
    buckling_failure_colors = ["olive", "seagreen", "maroon"]
    buckling_failure_labels = ["Fixed End Buckling", "Free End Buckling", "Web Buckling"]
    draw_other_BMDs(bm_tuples, bm_zeroes, buckling_failures, BRIDGE_LENGTH, buckling_failure_labels, buckling_failure_colors)
    plt.show()
    plt.clf()
    # Drawing the shear failures
    v_failures = V_fail(bridge, TauU, TauG)
    v_mat_fail = v_failures[0]
    v_top_fail = v_failures[1]
    shear_failures = [v_mat_fail, v_top_fail]
    shear_failure_labels = ["Material Shear", "Top Glue Shear"]
    shear_failure_colors = ["darkslategray", "orange"]
    if len(v_failures) > 2:
        v_bot_fail = v_failures[2]
        shear_failures.append(v_bot_fail)
        shear_failure_labels.append("Bottom Glue Shear")
        shear_failure_colors.append("aquamarine")  
    draw_other_SFDs(sf_tuples, shear_failures, BRIDGE_LENGTH, shear_failure_labels, shear_failure_colors)
    plt.show()
    plt.clf()
    v_buckling = VFailBuck(bridge, E, mu)
    shear_buck_failures = [v_buckling]
    s_buck_failure_labels = ["Shear Buckling"]
    s_buck_failure_colors = ["turquoise"]
    draw_other_SFDs(sf_tuples, shear_buck_failures, BRIDGE_LENGTH, s_buck_failure_labels, s_buck_failure_colors)
    plt.show()
    plt.clf()
    