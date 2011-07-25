## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
from __future__ import division

import numpy as np

6, 12, 24, 36, 48


vflux = {'2':(1/np.array([6, 12, 24, 36, 48]),
            np.array([(25.1077957032-2.47363376041j ),
                      (25.5014203462-1.32217653023j ),  
                      (25.7279169808-0.682546597138j), 
                      (25.8128815732-0.462930527586j), 
                      (25.8594020855-0.350742810606j)])),
         '2r':(1/np.array([6, 12, 24, 36, 48]),
               np.array([(25.437375959-2.47363376041j)  ,
                         (26.0343236965-1.32217653023j) ,  
                         (26.0267801945-0.682546597138j), 
                         (26.0226950991-0.462930527586j), 
                         (26.0212909762-0.350742810606j)])),

         '1':(1/np.array([6, 12, 24, 36, 48, 60, 72]),
            np.array([(27.2156108969-8.60337864309j ),  
                      (26.1174481478-4.69239069963j ),  
                      (25.6402182616-2.54217936348j ),  
                      (25.5353143779-1.75451501832j ),  
                      (25.6433969976-1.34895338793j ),
                      (25.8374635886-1.05092561784j ),
                      (25.809754997-0.888779820837j ),
                      ])),
         '1r':(1/np.array([6, 12, 24, 36, 48, 60, 72]),
               np.array([(29.7354509432-8.60337864309j),  
                         (27.598625396-4.69239069963j) ,  
                         (26.6933848464-2.54217936348j),  
                         (26.308038438-1.75451501832j) ,  
                         (26.2541326666-1.34895338793j),
                         (26.3222271753-1.05092561784j),
                         (26.2236911047-0.888779820837j),
                         ]))
            }

sflux = {'2':(1/np.array([6, 12, 24, 36, 48]),
            np.array([26.3586862457,
                      25.9114476977,  
                      25.9104357159, 
                      25.9325966765, 
                      25.9488038482])), 
         '1':(1/np.array([6, 12, 24, 36, 48]),
            np.array([32.1317821605,  
                      29.1373741302,  
                      27.3175968117,  
                      26.7039186893,  
                      26.5438681478]))
         }

voltage = {'2':(1/np.array([6, 12, 24, 36, 48]),
              np.array([(-26.8605732334+360.728686502j) ,
                        (-26.0024339562+866.66313718j)  ,  
                        (-25.9434293989+1162.24478969j) , 
                        (-25.9516416163+2137.09320662j) , 
                        (-25.9436917759+2836.97235168j) ])), 
           '1':(1/np.array([6, 12, 24, 36, 48]),
              np.array([(-31.5214674006+69.7869839465j),  
                        (-28.5605292573+205.076568969j),  
                        (-26.5259449392+391.725220001j),  
                        (-26.087188108+688.722412185j) ,  
                        (-26.029401834+972.238119062j) ]))
           }
