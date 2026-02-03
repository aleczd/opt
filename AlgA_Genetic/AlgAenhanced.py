############
############ ALTHOUGH I GIVE YOU THIS TEMPLATE PROGRAM WITH THE NAME 'skeleton.py', 
############ YOU CAN RENAME IT TO ANYTHING YOU LIKE. HOWEVER, FOR THE PURPOSES OF 
############ THE EXPLANATION IN THESE COMMENTS, I ASSUME THAT THIS PROGRAM IS STILL 
############ CALLED 'skeleton.py'.
############
############ IF YOU WISH TO IMPORT STANDARD MODULES, YOU CAN ADD THEM AFTER THOSE BELOW.
############ NOTE THAT YOU ARE NOT ALLOWED TO IMPORT ANY NON-STANDARD MODULES! TO SEE
############ THE STANDARD MODULES, TAKE A LOOK IN 'validate_before_handin.py'.
############
############ DO NOT INCLUDE ANY COMMENTS ON A LINE WHERE YOU IMPORT A MODULE.
############

import os
import sys
import time
import random
import math

############ START OF SECTOR 1 (IGNORE THIS COMMENT)
############
############ NOW PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############ BY 'DO NOT TOUCH' I REALLY MEAN THIS. EVEN CHANGING THE SYNTAX, BY
############ ADDING SPACES OR COMMENTS OR LINE RETURNS AND SO ON, COULD MEAN THAT
############ CODES MIGHT NOT RUN WHEN I RUN THEM!
############

def read_file_into_string(input_file, ord_range):
    the_file = open(input_file, 'r') 
    current_char = the_file.read(1) 
    file_string = ""
    length = len(ord_range)
    while current_char != "":
        i = 0
        while i < length:
            if ord(current_char) >= ord_range[i][0] and ord(current_char) <= ord_range[i][1]:
                file_string = file_string + current_char
                i = length
            else:
                i = i + 1
        current_char = the_file.read(1)
    the_file.close()
    return file_string

def remove_all_spaces(the_string):
    length = len(the_string)
    new_string = ""
    for i in range(length):
        if the_string[i] != " ":
            new_string = new_string + the_string[i]
    return new_string

def integerize(the_string):
    length = len(the_string)
    stripped_string = "0"
    for i in range(0, length):
        if ord(the_string[i]) >= 48 and ord(the_string[i]) <= 57:
            stripped_string = stripped_string + the_string[i]
    resulting_int = int(stripped_string)
    return resulting_int

def convert_to_list_of_int(the_string):
    list_of_integers = []
    location = 0
    finished = False
    while finished == False:
        found_comma = the_string.find(',', location)
        if found_comma == -1:
            finished = True
        else:
            list_of_integers.append(integerize(the_string[location:found_comma]))
            location = found_comma + 1
            if the_string[location:location + 5] == "NOTE=":
                finished = True
    return list_of_integers

def build_distance_matrix(num_cities, distances, city_format):
    dist_matrix = []
    i = 0
    if city_format == "full":
        for j in range(num_cities):
            row = []
            for k in range(0, num_cities):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    elif city_format == "upper_tri":
        for j in range(0, num_cities):
            row = []
            for k in range(j):
                row.append(0)
            for k in range(num_cities - j):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    else:
        for j in range(0, num_cities):
            row = []
            for k in range(j + 1):
                row.append(0)
            for k in range(0, num_cities - (j + 1)):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    if city_format == "upper_tri" or city_format == "strict_upper_tri":
        for i in range(0, num_cities):
            for j in range(0, num_cities):
                if i > j:
                    dist_matrix[i][j] = dist_matrix[j][i]
    return dist_matrix

def read_in_algorithm_codes_and_tariffs(alg_codes_file):
    flag = "good"
    code_dictionary = {}   
    tariff_dictionary = {}  
    if not os.path.exists(alg_codes_file):
        flag = "not_exist"  
        return code_dictionary, tariff_dictionary, flag
    ord_range = [[32, 126]]
    file_string = read_file_into_string(alg_codes_file, ord_range)  
    location = 0
    EOF = False
    list_of_items = []  
    while EOF == False: 
        found_comma = file_string.find(",", location)
        if found_comma == -1:
            EOF = True
            sandwich = file_string[location:]
        else:
            sandwich = file_string[location:found_comma]
            location = found_comma + 1
        list_of_items.append(sandwich)
    third_length = int(len(list_of_items)/3)
    for i in range(third_length):
        code_dictionary[list_of_items[3 * i]] = list_of_items[3 * i + 1]
        tariff_dictionary[list_of_items[3 * i]] = int(list_of_items[3 * i + 2])
    return code_dictionary, tariff_dictionary, flag

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ THE RESERVED VARIABLE 'input_file' IS THE CITY FILE UNDER CONSIDERATION.
############
############ IT CAN BE SUPPLIED BY SETTING THE VARIABLE BELOW OR VIA A COMMAND-LINE
############ EXECUTION OF THE FORM 'python skeleton.py city_file.txt'. WHEN SUPPLYING
############ THE CITY FILE VIA A COMMAND-LINE EXECUTION, ANY ASSIGNMENT OF THE VARIABLE
############ 'input_file' IN THE LINE BELOW iS SUPPRESSED.
############
############ IT IS ASSUMED THAT THIS PROGRAM 'skeleton.py' SITS IN A FOLDER THE NAME OF
############ WHICH IS YOUR USER-NAME, E.G., 'abcd12', WHICH IN TURN SITS IN ANOTHER
############ FOLDER. IN THIS OTHER FOLDER IS THE FOLDER 'city-files' AND NO MATTER HOW
############ THE NAME OF THE CITY FILE IS SUPPLIED TO THIS PROGRAM, IT IS ASSUMED THAT 
############ THE CITY FILE IS IN THE FOLDER 'city-files'.
############
############ END OF SECTOR 1 (IGNORE THIS COMMENT)

input_file = "AISearchfile012.txt"

############ START OF SECTOR 2 (IGNORE THIS COMMENT)
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if len(sys.argv) > 1:
    input_file = sys.argv[1]

############ END OF SECTOR 2 (IGNORE THIS COMMENT)
path_for_city_files = "../city-files"
############ START OF SECTOR 3 (IGNORE THIS COMMENT)
    
if os.path.isfile(path_for_city_files + "/" + input_file):
    ord_range = [[32, 126]]
    file_string = read_file_into_string(path_for_city_files + "/" + input_file, ord_range)
    file_string = remove_all_spaces(file_string)
    print("I have found and read the input file " + input_file + ":")
else:
    print("*** error: The city file " + input_file + " does not exist in the city-file folder.")
    sys.exit()

location = file_string.find("SIZE=")
if location == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
comma = file_string.find(",", location)
if comma == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
num_cities_as_string = file_string[location + 5:comma]
num_cities = integerize(num_cities_as_string)
print("   the number of cities is stored in 'num_cities' and is " + str(num_cities))

comma = comma + 1
stripped_file_string = file_string[comma:]
distances = convert_to_list_of_int(stripped_file_string)

counted_distances = len(distances)
if counted_distances == num_cities * num_cities:
    city_format = "full"
elif counted_distances == (num_cities * (num_cities + 1))/2:
    city_format = "upper_tri"
elif counted_distances == (num_cities * (num_cities - 1))/2:
    city_format = "strict_upper_tri"
else:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

dist_matrix = build_distance_matrix(num_cities, distances, city_format)
print("   the distance matrix 'dist_matrix' has been built.")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ YOU NOW HAVE THE NUMBER OF CITIES STORED IN THE INTEGER VARIABLE 'num_cities'
############ AND THE TWO_DIMENSIONAL MATRIX 'dist_matrix' HOLDS THE INTEGER CITY-TO-CITY 
############ DISTANCES SO THAT 'dist_matrix[i][j]' IS THE DISTANCE FROM CITY 'i' TO CITY 'j'.
############ BOTH 'num_cities' AND 'dist_matrix' ARE RESERVED VARIABLES AND SHOULD FEED
############ INTO YOUR IMPLEMENTATIONS.
############
############ THERE NOW FOLLOWS CODE THAT READS THE ALGORITHM CODES AND TARIFFS FROM
############ THE TEXT-FILE 'alg_codes_and_tariffs.txt' INTO THE RESERVED DICTIONARIES
############ 'code_dictionary' AND 'tariff_dictionary'. DO NOT AMEND THIS CODE!
############ THE TEXT FILE 'alg_codes_and_tariffs.txt' SHOULD BE IN THE SAME FOLDER AS
############ THE FOLDER 'city-files' AND THE FOLDER WHOSE NAME IS YOUR USER-NAME.
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############
############ END OF SECTOR 3 (IGNORE THIS COMMENT)

############ START OF SECTOR 4 (IGNORE THIS COMMENT)
path_for_alg_codes_and_tariffs = "../alg_codes_and_tariffs.txt"
############ END OF SECTOR 4 (IGNORE THIS COMMENT)

############ START OF SECTOR 5 (IGNORE THIS COMMENT)
code_dictionary, tariff_dictionary, flag = read_in_algorithm_codes_and_tariffs(path_for_alg_codes_and_tariffs)

if flag != "good":
    print("*** error: The text file 'alg_codes_and_tariffs.txt' does not exist.")
    sys.exit()

print("The codes and tariffs have been read from 'alg_codes_and_tariffs.txt':")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY! SORRY TO GO ON ABOUT THIS BUT YOU NEED TO BE 
############ AWARE OF THIS FACT!
############
############ YOU NOW NEED TO SUPPLY SOME PARAMETERS.
############
############ THE RESERVED STRING VARIABLE 'my_user_name' SHOULD BE SET AT YOUR
############ USER-NAME, E.G., "abcd12"
############
############ END OF SECTOR 5 (IGNORE THIS COMMENT)

my_user_name = "pwtn63"

############ START OF SECTOR 6 (IGNORE THIS COMMENT)
############
############ YOU CAN SUPPLY, IF YOU WANT, YOUR FULL NAME. THIS IS NOT USED AT ALL BUT SERVES AS
############ AN EXTRA CHECK THAT THIS FILE BELONGS TO YOU. IF YOU DO NOT WANT TO SUPPLY YOUR
############ NAME THEN EITHER SET THE STRING VARIABLES 'my_first_name' AND 'my_last_name' AT 
############ SOMETHING LIKE "Mickey" AND "Mouse" OR AS THE EMPTY STRING (AS THEY ARE NOW;
############ BUT PLEASE ENSURE THAT THE RESERVED VARIABLES 'my_first_name' AND 'my_last_name'
############ ARE SET AT SOMETHING).
############
############ END OF SECTOR 6 (IGNORE THIS COMMENT)

my_first_name = "Alec"
my_last_name = "Durgheu"

############ START OF SECTOR 7 (IGNORE THIS COMMENT)
############
############ YOU NEED TO SUPPLY THE ALGORITHM CODE IN THE RESERVED STRING VARIABLE 'algorithm_code'
############ FOR THE ALGORITHM YOU ARE IMPLEMENTING. IT NEEDS TO BE A LEGAL CODE FROM THE TEXT-FILE
############ 'alg_codes_and_tariffs.txt' (READ THIS FILE TO SEE THE CODES).
############
############ END OF SECTOR 7 (IGNORE THIS COMMENT)

algorithm_code = "GA"

############ START OF SECTOR 8 (IGNORE THIS COMMENT)
############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS STARTING
############ 'HAVE YOU TOUCHED ...'
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if not algorithm_code in code_dictionary:
    print("*** error: the algorithm code " + algorithm_code + " is illegal")
    sys.exit()
print("   your algorithm code is legal and is " + algorithm_code + " -" + code_dictionary[algorithm_code] + ".")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY! SORRY TO GO ON ABOUT THIS BUT YOU NEED TO BE 
############ AWARE OF THIS FACT!
############
############ YOU CAN ADD A NOTE THAT WILL BE ADDED AT THE END OF THE RESULTING TOUR FILE IF YOU LIKE,
############ E.G., "in my basic greedy search, I broke ties by always visiting the first 
############ city found" BY USING THE RESERVED STRING VARIABLE 'added_note' OR LEAVE IT EMPTY
############ IF YOU WISH. THIS HAS NO EFFECT ON MARKS BUT HELPS YOU TO REMEMBER THINGS ABOUT
############ YOUR TOUR THAT YOU MIGHT BE INTERESTED IN LATER.
############
############ END OF SECTOR 8 (IGNORE THIS COMMENT)

added_note = ""

############
############ NOW YOUR CODE SHOULD BEGIN.
############


########################################################################################################################################

# Given Variables
# num_cities = The number of cities in the instance of TSP
# dist_matrix = dist_matrix[i][j] gives the iinteger distance from city i to j, recall the graph is a clique of size V


"""
------------------------------------------- HYPERPARAMETERS -------------------------------------------------------
"""



max_it = 99999999999999999999              # the maximum number of iterations (arbitrary due to use of a killswitch)
population_size = 10                       # the size of the population
prob_mutation = 0.3                        # the probability of a mutation
threshold_time = 55                        # controls the duration of the search (in seconds) (for the killswitch)

elite_size = 4                             # the size of the elite group from the population
harsh_mute = 0.1                           # the probability of a mutation being a harsh mutation



"""
---------------------------------------------- GLOBAL VARIABLES ----------------------------------------------
"""


P = [[]] * population_size                 # stores the population of individuals at time t
P_fitness = [[]] * population_size         # stores the fitness of each individual at time t
newP = []                                  # used to calculate the individuals at time t+1
newP_fitness = []                          # stores the fitness of the individuals at time t+1
init_time = time.time()                    # starts the timer for the killswitch routine
tau = 0                                    # constant used in the fitness function
b_tour = []                                # contains the best tour found so far
b_tour_length = 0                          # contains the length of the best tour found so far

elites = []                                # stores the elite group from the population at time t
elites_fitness = []                        # stores the fitness of the elite group                       
timer_t = 0                                # tracks the in-script time of important events, eg checkpoints
devol = 2                             # the intensity of the devolve function
cook_t = 1000                              # the in-script time allowed before identifying as stuck in a local optimum
iters = 0                                  # keeps track of the total number of useful iterations (in 'in-script' time)
saving = False                             # a flag used to ignore certain in-script time iterations


"""
--------------------------------------------- FUNCTIONS-----------------------------------------------
"""





"""
- This function returns random tours and is predominantly used to initialise the algorithm
"""
def random_tour():
    random_tour = list(range(0,num_cities))
    random.shuffle(random_tour)
    return random_tour








"""
- Given an individual, this function will return the tour length
"""
def calculate_length(ind):
    length = 0
    for i in range(0,num_cities-1):
        length += dist_matrix[ind[i]][ind[i+1]]
    return length + dist_matrix[ind[-1]][ind[0]]








"""
- This function is responsible for calculating the number tau for use in the fitness function
- It is used to set the global value of tau in initialise_population()
"""
def get_tau():
    tau = 0
    for i in range(0, len(dist_matrix)):
        tau += math.floor(1.5 * (sum(dist_matrix[i])/num_cities))
    return tau









"""
- This function calculates the fitness of a given individual using tau
- It will ensure that the minimum fitness value is (lower) bounded by 1 in the case
    of a particularly ineffective tour 
- It returns this fitness value
"""
def fitness(ind):
    trial = tau - calculate_length(ind)
    if trial < 1:
        return 1         
    return trial










"""
- This function calculates the best tour from the population
- It returns the best tour, as well as the length of the best tour
"""
def best_tour():
    best = 0
    for i in range(1, population_size):
        if P_fitness[i] > P_fitness[best]:
            best = i
    return (P[best], calculate_length(P[best]))











"""
- This is the procedure that initialises the genetic algorithm
- It randomly creates a population and determines its fitness, as well as the best tour so far
"""
def initialise_population():

    # Global Variables Editted
    global P
    global P_fitness
    global tau
    global b_tour
    global b_tour_length

    # Set Tau globally
    tau = get_tau()

    # Build the initial population
    for i in range(0, population_size):
        new_tour = random_tour()
        P[i] = new_tour
        P_fitness[i] = fitness(new_tour)

    # Find and set the best tour from the initial population
    b_tour, b_tour_length = best_tour()















"""
- This is essentially the roulette wheel function, using the random.choices method
- It selects from the elite group (elites), weighted by the fitness of the elites (elites_fitness)
- It returns the selected individual
"""
def select_individual():
    return random.choices(elites, weights = elites_fitness, k = 1)[0] 










"""
- This function takes a given encoding of a tour and mutates it with a certain chance (prob_mutation)
- The mutation protocol:
    - if a harsh mutation occurs, it will swap two elements (as in the basic algorithm)
    - otherwise, the regular mutation would find two breakpoints and reverse the tour between these breakpoints
- We return the mutated individual if a mutation occurs, otherwise we return the input chromosome (no mutation)
"""
def mutate(x):
    if random.random() < prob_mutation:
        if random.random() < harsh_mute:

            # Swap two elements randomly (harsh mutation - 4 changes)
            i1 = random.randint(0, num_cities-1)
            i2 = random.randint(0, num_cities-1)
            while i2 == i1:
                i2 = random.randint(0, num_cities-1)
            temp = x[i1]
            x[i1] = x[i2]
            x[i2] = temp
        else:

            # Middle reverse (general mutation - 2 changes)
            start = random.randint(1, num_cities-3)
            end = random.randint(start + 1, num_cities-2)
            x = x[:start] + x[end:start-1:-1] + x[end+1:]

    return x







"""
- This function handles the crossover of two parent tours (x, y)
- It uses the crossover technique presented in the lectures, a one point crossover
- The two children are potentially mutated before returning the one with the greater fitness
- More specifically, we return the fitter child tour and its fitness
"""
def reproduce(x, y):

    # Randomly slice the genes
    slice_index = random.randint(1, num_cities-1)
    x_suff = x[slice_index:]
    x_pref = x[:slice_index]
    y_suff = y[slice_index:]
    y_pref = y[:slice_index]

    # Ammend children
    y_key = 0
    x_key = 0
    for i in range(0, len(y_suff)):
        while y_suff[i] in x_pref:
            y_suff[i] = y_pref[y_key]
            y_key += 1
        while x_suff[i] in y_pref:
            x_suff[i] = x_pref[x_key]
            x_key += 1

    # Apply mutation to children
    z_x = mutate(x_pref + y_suff)
    z_y = mutate(y_pref + x_suff)

    # Determine which child will be returned
    z_x_fitness = fitness(z_x)
    z_y_fitness = fitness(z_y)
    if z_x_fitness >= z_y_fitness:
        return (z_x, z_x_fitness)
    return (z_y, z_y_fitness)







"""
- This function keeps track of the forceful termination condition (killswitch)
- It returns True if the program should terminate, as well as the timestamp value
"""
def timeout():
    val = time.time() - init_time
    if val > threshold_time:
        return (True, val)
    return (False, val)









"""
- This function will determine the elite group from a given population (copy_P) and
    its fitness (copy_fitness)
- It edits the global variables to store the elite group at time t
"""
def find_elites(copy_fitness, copy_P):

    # Global Variables Editted
    global elites
    global elites_fitness

    # Reset Globals
    elites = []
    elites_fitness = []

    # Repeatedly select the chromosomes with the highest fitness values
    for i in range(0, elite_size):
        max_key = 0
        for j in range(1, len(copy_fitness)):
            if copy_fitness[j] > copy_fitness[max_key]:
                max_key = j
        elites.append(copy_P[max_key])
        elites_fitness.append(copy_fitness[max_key])
        copy_fitness.pop(max_key)
        copy_P.pop(max_key)















"""
- This function takes a given tour and devolves it according 
    to the intensity (devol)
- We essentially perform a "left shift" on the indeces selected to swap
    - this means the first position ends up in the last, 
    - the last position ends up in the penultimate position
    - etc...
    - devol specifies the number of indeces, and hence the number of
        swaps performed
- We return the devolved chromosome
"""
def devolve(x):
    indeces = random.sample(list(range(0,num_cities)), k = devol)
    for i in range(0, len(indeces)-1):
        temp = x[indeces[i]]
        x[indeces[i]] = x[indeces[i+1]]
        x[indeces[i+1]] = temp
    return x










"""
- This function recreates an ideal scenario for the genetic algorithm to 
    retry searching
- It uses the best tour found so far and devolves it to rebuild a new population
- The devol increases cyclically to allow for more branches to be reconsidered
- The algorithm can 'cook' (search) for an adjusted duration of time as the search 
    gets more sophisticated (by using iter)
- It returns the new population and its fitness for the algorithm to proceed
"""
def checkpoint():

    # Global Variables Editted
    global devol
    global cook_t

    # Essentially, make everything to be the best tour we found but 
    # "mutate somewhat harshly" (devolve)
    cP = []
    cP_fitness = []

    # Rebuild checkpoint population
    for i in range(0, population_size):

        # Create copy of best tour
        cb_tour = []
        for i in b_tour:
            cb_tour.append(i)

        # Devolve
        new = devolve(cb_tour)
        cP.append(new)
        cP_fitness.append(fitness(new))

    # Increase devol to be at most half of the number of cities
    if devol < num_cities//2:
        devol += 1
    else:
        devol = 2
    
    # Allow the search to continue for longer depending on our current depth
    cook_t = iters + 1000

    return (cP, cP_fitness)












"""
- This is the main procedure for the Genetic algorithm
- It uses all of the funtions defined above and executes the algorithm
- It takes the maximum number of iterations as a parameter (max_it - though we have global access)
"""
def GA(max_it):

    # Global Variables Editted
    global P
    global P_fitness
    global newP
    global newP_fitness
    global b_tour
    global b_tour_length

    global timer_t
    global devol
    global cook_t
    global iters
    global saving

    # Initialise the algorithm
    initialise_population()

    # The main loop
    t = 0
    while t < max_it:
        newP = []
        newP_fitness = []
        
        # Find the elite group before reproducing the next generation
        find_elites(P_fitness, P)

        for i in range(0, population_size):

            # Select two parents, weighted by their fitness
            X = select_individual()
            Y = select_individual()

            # Obtain the child from the crossover and mutation operations
            Z, Z_fitness = reproduce(X, Y)

            # Add this child to the next generation population
            newP.append(Z)
            newP_fitness.append(Z_fitness)
        
        # If there has been no new discovery for a while, run checkpoint()
        # The saving flag indicates that we are current checkpointing and
        #   iters should not be updated
        if (t - timer_t > cook_t):
            newP, newP_fitness = checkpoint()
            timer_t = t
            saving = True
    
        # After the iteration at time t...
        # Swap the data in preparation for time t+1
        temp = newP
        P = temp

        temp2 = newP_fitness
        P_fitness = temp2

        # Calculate the best tour from this population and update the global best tour thus far if need be
        newBest, newBestLength = best_tour()
        if newBestLength < b_tour_length:
            b_tour = newBest
            b_tour_length = newBestLength

            # Check the saving flag to see if we should update iters
            if not saving:
                iters += t - timer_t
            
            # Reset the various chekpointing methods to their default
            timer_t = t
            devol = 2
            cook_t = 1000
            saving = False

            #print("New best overall: ", b_tour_length)

        # Increase the time t to t+1 and check if the program should forcefully terminate
        t += 1
        occ, val = timeout()
        if occ:
            break

    # Return the best tour found so far (Note: it is also accessible globally)
    return (b_tour, b_tour_length)





"""
------------------------------------------------ The Script ----------------------------------------------------
"""



"""
Run the algorithm and assign the variables correctly
"""
GA(max_it)
tour = b_tour
tour_length = b_tour_length




"""
-------------------------------------------------------------------------------------------------------------------------
"""









############ START OF SECTOR 9 (IGNORE THIS COMMENT)
############
############ YOUR CODE SHOULD NOW BE COMPLETE AND WHEN EXECUTION OF THIS PROGRAM 'skeleton.py'
############ REACHES THIS POINT, YOU SHOULD HAVE COMPUTED A TOUR IN THE RESERVED LIST VARIABLE 'tour', 
############ WHICH HOLDS A LIST OF THE INTEGERS FROM {0, 1, ..., 'num_cities' - 1} SO THAT EVERY INTEGER
############ APPEARS EXACTLY ONCE, AND YOU SHOULD ALSO HOLD THE LENGTH OF THIS TOUR IN THE RESERVED
############ INTEGER VARIABLE 'tour_length'.
############
############ YOUR TOUR WILL BE PACKAGED IN A TOUR FILE OF THE APPROPRIATE FORMAT AND THIS TOUR FILE'S,
############ NAME WILL BE A MIX OF THE NAME OF THE CITY FILE, THE NAME OF THIS PROGRAM AND THE
############ CURRENT DATE AND TIME. SO, EVERY SUCCESSFUL EXECUTION GIVES A TOUR FILE WITH A UNIQUE
############ NAME AND YOU CAN RENAME THE ONES YOU WANT TO KEEP LATER.
############
############ DO NOT TOUCH OR ALTER THE CODE BELOW THIS POINT! YOU HAVE BEEN WARNED!
############

flag = "good"
length = len(tour)
for i in range(0, length):
    if isinstance(tour[i], int) == False:
        flag = "bad"
    else:
        tour[i] = int(tour[i])
if flag == "bad":
    print("*** error: Your tour contains non-integer values.")
    sys.exit()
if isinstance(tour_length, int) == False:
    print("*** error: The tour-length is a non-integer value.")
    sys.exit()
tour_length = int(tour_length)
if len(tour) != num_cities:
    print("*** error: The tour does not consist of " + str(num_cities) + " cities as there are, in fact, " + str(len(tour)) + ".")
    sys.exit()
flag = "good"
for i in range(0, num_cities):
    if not i in tour:
        flag = "bad"
if flag == "bad":
    print("*** error: Your tour has illegal or repeated city names.")
    sys.exit()
check_tour_length = 0
for i in range(0, num_cities - 1):
    check_tour_length = check_tour_length + dist_matrix[tour[i]][tour[i + 1]]
check_tour_length = check_tour_length + dist_matrix[tour[num_cities - 1]][tour[0]]
if tour_length != check_tour_length:
    flag = print("*** error: The length of your tour is not " + str(tour_length) + "; it is actually " + str(check_tour_length) + ".")
    sys.exit()
print("You, user " + my_user_name + ", have successfully built a tour of length " + str(tour_length) + "!")

local_time = time.asctime(time.localtime(time.time()))
output_file_time = local_time[4:7] + local_time[8:10] + local_time[11:13] + local_time[14:16] + local_time[17:19]
output_file_time = output_file_time.replace(" ", "0")
script_name = os.path.basename(sys.argv[0])
if len(sys.argv) > 2:
    output_file_time = sys.argv[2]
output_file_name = script_name[0:len(script_name) - 3] + "_" + input_file[0:len(input_file) - 4] + "_" + output_file_time + ".txt"

f = open(output_file_name,'w')
f.write("USER = " + my_user_name + " (" + my_first_name + " " + my_last_name + "),\n")
f.write("ALGORITHM CODE = " + algorithm_code + ", NAME OF CITY-FILE = " + input_file + ",\n")
f.write("SIZE = " + str(num_cities) + ", TOUR LENGTH = " + str(tour_length) + ",\n")
f.write(str(tour[0]))
for i in range(1,num_cities):
    f.write("," + str(tour[i]))
f.write(",\nNOTE = " + added_note)
f.close()
print("I have successfully written your tour to the tour file:\n   " + output_file_name + ".")

############ END OF SECTOR 9 (IGNORE THIS COMMENT)
    
    











    


