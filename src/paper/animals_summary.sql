-- Animals used in behavior
SELECT animalid,sex,construct_name,dob,ephysid FROM mice.animals WHERE animalid IN 
('19-MI16506-02',
'19-MI16506-03',
'19-MI17966-03',
'19-MI17966-04',
'20-MI10996-01',
'20-MI10996-02');


-- Animals used in electrophysiology with units with baseline+salicylate conditions
SELECT animalid,sex,construct_name,dob,ephysid FROM mice.animals WHERE ephysid IN 
('M25', 'M26', 'M27', 'M29', 
 'M30', 'M42', 'M43', 
'M71',  'M72', 'M74');

-- Counting the sexes from the data above
SELECT sex, count(animalid) FROM 
(SELECT animalid,sex FROM mice.animals WHERE ephysid IN 
('M25', 'M26', 'M27', 'M29', 
 'M30', 'M42', 'M43', 
'M71',  'M72', 'M74')
 ) AS tb
 GROUP BY sex;
 
 -- Animals used in both behavior and in electrophysiology with units with baseline+salicylate conditions
SELECT animalid,sex,construct_name,dob,ephysid FROM mice.animals WHERE animalid IN 
('19-MI16506-02',
'19-MI16506-03',
'19-MI17966-03',
'19-MI17966-04',
'20-MI10996-01',
'20-MI10996-02')
AND 
ephysid IN 
('M25', 'M26', 'M27', 'M29', 
 'M30', 'M42', 'M43', 
'M71',  'M72', 'M74');

-- [NOT IN DATABASE!!] animals used in imaging 
SELECT animalid,sex,construct_name,dob,ephysid FROM mice.animals WHERE animalid IN 
('15-MI13723-01',-- [NOT IN DATABASE!!] 
'15-MI13723-02');-- [NOT IN DATABASE!!] 
