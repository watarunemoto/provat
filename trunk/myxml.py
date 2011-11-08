
from procrun import fileOpen
import re

class myXmlParser :
	def __init__(self): self.handler = None
	def setContentHandler(self,handler) : self.handler = handler
	def parse(self,filename) :
		f = fileOpen(filename,'r')
		alltext, tokens, istag = "", [], []
		for l in f.readlines():
			if re.compile("^$").search(l) : continue
			alltext = alltext + l
		index = 0
		while 1:
			m = re.compile("<[^<]*>").search(alltext[index:])
			if not m : break
			#print alltext[index:index+3], m.group()
			if m.start() > 0 :
				tokens.append(alltext[index:index+m.start()])
				istag.append(0)
			tokens.append( alltext[index+m.start()+1 : index+m.end()-1] )
			istag.append(1)
			index = index + m.end()
		#for t in tokens : print "TOKEN",t
		for ti in range(len(tokens)) :
			if istag[ti] == 1 :
				if tokens[ti][0] == '/' : self.handler.endElement(tokens[ti][1:])
				else : self.handler.startElement(tokens[ti],"")
			else : self.handler.characters(tokens[ti])
		self.handler.endDocument()
