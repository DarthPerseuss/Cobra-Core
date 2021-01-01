# Mining RNA epigenetics data from the web

import os
from selenium import webdriver
import pandas as pd
import numpy as np
import time



nonCore_df = pd.read_excel("genes.xlsx", header=None)
core_df = pd.read_excel("Core Genes.xlsx", header=None)
for i in core_df.values:
    nonCore_df = nonCore_df[nonCore_df != i]
nonCore_df = nonCore_df.dropna()


# get the path of ChromeDriverServer
direct = os.path.dirname(__file__)
chrome_driver_path = direct + "\chromedriver.exe"

# create a new Chrome session
driver = webdriver.Chrome(chrome_driver_path)
driver.implicitly_wait(5)
driver.maximize_window()

# navigate to the application home page
driver.get("http://lulab.life.tsinghua.edu.cn/clipdb/targetsearch.php")

# selects yeast from the dropbox
species_field = driver.find_element_by_id("select-trigger-genes")
species_select = driver.find_element_by_id("select-content-genes")
binding_targets = driver.find_element_by_xpath('//*[@id="success-table-genes"]/tbody')
hint = driver.find_element_by_id("genesHint")

# initialize the search protocol
species_field.click()
species_select.find_element_by_id("3").click()
search_field = driver.find_element_by_id("inputgenesmanul")
search_field.clear()
search_field.send_keys("YAL012W")
driver.find_element_by_id("search-genes").click()

for i in binding_targets.find_elements_by_xpath('.//td'):
    print(i.get_attribute('innerHTML'))

for i in range(8, len(nonCore_df)):
    search_field = driver.find_element_by_id("inputgenesmanul")  # input gene
    if i != len(nonCore_df):
        search_field.clear()
    else:
        pass
    search_field.send_keys((nonCore_df.iat[i, 0]))
    driver.find_element_by_id("search-genes").click()
    time.sleep(5)
    if binding_targets.is_displayed() is True:
        print(nonCore_df.iat[i, 0])
        for j in binding_targets.find_elements_by_xpath('.//td'):
            print(j.get_attribute('innerHTML'))

    else:
        print(nonCore_df.iat[i, 0] + " NaN")
