import csv
import pandas as pd

df = pd.read_csv('titanic-data.csv')

#Print basic statistics about the Titanic dataset
print(df.describe())

#Looking for youngest passenger
print(df.loc[df['Age'] == 1])

#Looking for oldest passenger
print(df.loc[df['Age'] == 92])

#Creating new cabin level column to tidy cabin column
def first_letter(string):
	return string[:1]

df['Cabin Letters'] = df['Cabin'].apply(first_letter)

#Grouping function of passenger age groups depending on whether they survived or not

def titanic_age(column1, survived):
	for i in range(0, 100, 10):
		print(str(i) + ' -', str(i+10), ': ', (len(df['PassengerId'][(df['Survived'] == survived) 
			& ((df[column1] >= i) & (df[column1] < i+10))])))

#Grouping function of passenger sex groups depending on whether they survived or not

def titanic_sex(column1, survived, sex):
	print(len(df['PassengerId'][(df['Survived'] == survived) & (df[column1] == sex)]))


#Grouping function of passenger class groups depending on whether they survived or not

def titanic_class(column1, survived):
	for i in range(1, 4):
		print(len(df['PassengerId'][(df['Survived'] == survived) & (df[column1] == i)]))



#Age groups divided by decades in relation to survivorship
print('\nSurvived, Age \n')
titanic_age('Age', 1)
print('\nDied, Age \n')
titanic_age('Age', 0)

# #Sex groups divided by decades in relation to survivorship
# print('\nSurvived, Sex \n')
# print('Male: ')
# titanic_sex('Sex', 1, 'male')
# print('Female: ')
# titanic_sex('Sex', 1, 'female')
# print('\nDied, Sex \n')
# print('Male: ')
# titanic_sex('Sex', 0, 'male')
# print('Female: ')
# titanic_sex('Sex', 0, 'female')

# #Class groups divided by decades in relation to survivorship
# print('\nSurvived, Class \n')
# pclass_survived = df.groupby('Pclass').sum()['Survived']
# print('\nDied, Class \n')
# titanic_class('Pclass', 0)

total_male_female = df['Sex'].value_counts()
print(total_male_female)

class_survival = df.groupby('Pclass').sum()['Survived']
print(class_survival)







