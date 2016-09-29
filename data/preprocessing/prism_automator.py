from selenium import webdriver

browser = webdriver.Firefox()

browser.get('http://tomcat:tomcat@localhost:8080/manager/html')

warfield = browser.find_element_by_name("deployWar")

button = browser.find_element_by_xpath("//input[@type='submit'][@value='Deploy']")

browser.get('http://tomcat:tomcat@localhost:8080/prism/#!/prism')

genome_file = browser.find_element_by_class_name("textfw")
