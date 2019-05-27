#!/usr/local/bin/python

import os
import sys 
import pickle
import os.path
import re
import csv
import urllib
from urllib.request import urlopen
from bs4 import BeautifulSoup
import pdb

# Reference websites
stks_website = "http://www.sse.com.cn/market/sseindex/indexlist/s/i000016/const_list.shtml"

stks_reference = "http://vip.stock.finance.sina.com.cn/corp/go.php/vCB_AllBulletin/stockid/" # then add str(stkid) and ".phtml"

# From here, we will first get a list of tags, from which to extract stock id's. This will be website_lists.

with urllib.request.urlopen(stks_website) as response:
   html = response.read()

website_soup = BeautifulSoup(html, "lxml") 

website_lists = website_soup.find_all('td', {'class': 'table3'})

website_lists = website_lists[:-1]

website_ids = []

for website_bstag in website_lists:
    
    website_string = str(website_bstag)

    website_start = website_string.find("COMPANY_CODE") + 13
    
    website_end = website_start + 6
    
    website_stkid = website_string[ website_start: website_end ]
    
    website_ids.append(website_stkid)


# ----------------------------------------------------------------------------------
# After getting the list of stock id's, we are for good done with the first website.

# Let's test it on one of them.
    
with urllib.request.urlopen(stks_reference + "603993" + ".phtml") as response:
   html = response.read()

your_soup = BeautifulSoup(html, "lxml") # good practice_under_score

# Use help() to check documentation.

your_lists = your_soup.find('div', {'class': 'datelist'}).find_all("a", {"target": "_blank"}) 

# Experimenting with the link of the first announcement:
your_string = str( your_lists[0] )

your_start = your_string.find("_blank") + 8

your_end = your_string.find("</a>")

your_title = your_string[your_start : your_end] # There, we have the title of the announcement.

# Let's get the date. First we need to get together the weblink for each individual announcement.

stks_announce = stks_reference[:stks_reference.find("cn") + 2] # http://vip.stock.finance.sina.com.cn

# Getting the second part of the website of the announcement from your_string. We also need to eliminate amp;

announce_start = your_string.find("/")

announce_mid = your_string.find("amp")

announce_end = your_string.find("target") - 2

announce_part = your_string[announce_start : announce_mid] + your_string[announce_mid + 4 : announce_end]

announce_overall = stks_announce + announce_part

with urllib.request.urlopen(announce_overall) as response:
   html = response.read()

his_soup = BeautifulSoup(html, 'lxml')

his_lists = his_soup.find_all("font", {"size" : "2"})

his_string = str(his_lists[0])

his_start = his_string.find(".PDF") - 18

his_end = his_string.find(".PDF") - 8

announce_date = his_string[his_start : his_end]

# Now, we need to see how to access the text and count the number of words. We won't retain the ones with strictly less than 50.

announce_content = ""

text_lists = his_soup.find_all("p")

for i in range( 0, len(text_lists) ):
    
    content = text_lists[i].text
    
    announce_content += content


# The loop will repeat exactly what we did up there for every stock. Don't run this yet. This will take ages.
for identity in website_ids:
    
    current_link = stks_reference + identity + ".phtml"
    
    with urllib.request.urlopen(current_link) as response:
       html = response.read()
    
    my_soup = BeautifulSoup(html, "lxml")

# Then put them back into the webpage.

# Then urllib and beautiful soup all of them.

# csv: stockid, article_number, date, headline

# find new pattern