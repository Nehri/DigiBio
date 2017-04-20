'''
Usage:
source venv/bin/activate
python automatic_wrappa_user.python
'''

import selenium
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
import time
import os
import argparse
import re
import urllib

pdb_files = ["1dm1", "1lhs"]

def main():
	for file in pdb_files:

		# Starts the webdriver
		driver = webdriver.Firefox()
		driver.get("http://www.wrappa.org/wrappa01/")
		src = driver.page_source

		# Moves past the terms and agreements page
		driver.find_element_by_name("termsAccepted").click()
		agree_button_xpath = '//input[@type="submit" and @value="agree"]'
		agree_element = driver.find_element_by_xpath(agree_button_xpath)
		agree_element.click()

		# Opens the wrappa analyzer
		driver.get("http://www.wrappa.org/wrappa01/wrappa")

		# Chooses the file to be analyzed
		choose_button_xpath = '//input[@type="file" and @name="pdbFileName"]'
		choose_element = driver.find_element_by_xpath(choose_button_xpath)
		choose_element.send_keys(os.path.dirname(os.path.abspath(file+".pdb"))+"/"+file+".pdb")

		# Uploads the pdb file
		upload_button_xpath = '//input[@type="submit" and @value="Upload"]'
		upload_element = driver.find_element_by_xpath(upload_button_xpath)
		upload_element.click()

		# Configures the default configuration
		configure_button_xpath = '//input[@type="submit" and @value="configure"]'
		configure_element = driver.find_element_by_xpath(configure_button_xpath)
		configure_element.click()

		# Analyzes the file
		analyze_button_xpath = '//input[@type="submit" and @value="analyze"]'
		analyze_element = driver.find_element_by_xpath(analyze_button_xpath)
		analyze_element.click()

		# Saves the output files
		bonds_link = driver.find_element_by_link_text("Bonds")
		bonds_link.click()
		bonds_src = driver.page_source
		with open(file+"_bonds_default.txt", "w+") as f:
			f.write(bonds_src)
		driver.back()

		wrappers_link = driver.find_element_by_link_text("Wrappers")
		wrappers_link.click()
		wrappers_src = driver.page_source
		with open(file+"_wrappers_default.txt", "w+") as f:
			f.write(wrappers_src)
		driver.back()

		anomaly_link = driver.find_element_by_link_text("Anomalies")
		anomaly_link.click()
		anomaly_src = driver.page_source
		with open(file+"_anomalies_default.txt", "w+") as f:
			f.write(anomaly_src)

		# Closes the browser window
		driver.close()
		
		# Ends the current firefox session
		driver.quit

if __name__ == '__main__':
	main()

