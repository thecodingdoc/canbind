CXXFLAGS := -O3 -Wall 
INCLUDES := -Iboost

all: clusterBindingSites weightByContacts calculatePValues

clusterBindingSites: clusterBindingSites.C
	g++ $(INCLUDES) $(CXXFLAGS) -o clusterBindingSites clusterBindingSites.C utils.C

weightByContacts: weightByContacts.C
	g++ $(INCLUDES) $(CXXFLAGS) -o weightByContacts weightByContacts.C utils.C

calculatePValues: calculatePValues.C
	g++ $(INCLUDES) $(CXXFLAGS) -o calculatePValues calculatePValues.C utils.C

clean:
	rm clusterBindingSites weightByContacts calculatePValues
