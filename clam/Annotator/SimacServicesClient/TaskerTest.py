#!/usr/bin/python

# Copyright (c) 2001-2006 MUSIC TECHNOLOGY GROUP (MTG)
#                         UNIVERSITAT POMPEU FABRA
#
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


from Tasker import Tasker, TaskerError 
import unittest
from cStringIO import StringIO
import os

taskfile=StringIO("""<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<Task>
  <Description>This is a really nice task</Description>
  <ContentLocator>http://localhost/SimacServices/ContentLocator</ContentLocator>
  <MetadataProvider>http://localhost/SimacServices/MetadataProvider</MetadataProvider>
  <IDList>
    <ID>69</ID>
  </IDList>
  <Descriptors>
    <Descriptor>Song::Artist</Descriptor>
    <Descriptor>Song::Title</Descriptor>
  </Descriptors>
</Task>""")

wrongtaskfile=StringIO("""<Task/>""")

def nullPrint(message):
	pass

class TaskerTest ( unittest.TestCase ) :
	def setUp(self):
		self.tasker = Tasker(nullPrint)

	def testRetrieveTask_nonExistingorIncorrectTaskFile(self):
		try:
			self.tasker.retrieveTask( "nonexisting.file" )
			self.fail( "A TaskerError should have been thrown given a non-existing task file" )
		except TaskerError:
			pass
		try:
			self.tasker.retrieveTask( StringIO("<I am a wrong XML file!!!") )
			self.fail( "A TaskerError should have been thrown given a non-existing task file" )
		except TaskerError:
			pass

	def testExtractParameters_correctFile(self):
		self.tasker.setParameters( taskfile, "Project" )
		ids, descs, cl, mp, description = self.tasker.extractParameters() 

		self.assertEquals( ids, ['69'] )
		self.assertEquals( descs, ['Song::Artist','Song::Title'] )
		self.assertEquals( cl,'http://localhost/SimacServices/ContentLocator' )
		self.assertEquals( mp,'http://localhost/SimacServices/MetadataProvider' )
		self.assertEquals( description,'This is a really nice task' )

	def testExtractParameters_wrongFile(self):
		self.tasker.setParameters( wrongtaskfile, "Project" )
		try:
			ids, descs, cl, mp = self.tasker.extractParameters() 
		except TaskerError:
			pass
		else:
			self.fail("A TaskError should have been thrown")
	
	def testCreateFile(self):
		self.tasker.createFile('erase', '.me', 'this is a content')
		file = open('erase.me')
		self.assertEquals( file.read(), 'this is a content')
		os.remove('erase.me')

	def testDownloadSong_noLocationsSpecified(self):
		result=self.tasker.downloadSong("")
		self.assertEquals( result, None)

	def testDownloadSong_noValidLocations(self):
		result=self.tasker.downloadSong( "invalidpath\noh://my.god" )
		self.assertEquals( result, None)

if __name__ == "__main__" :
	suite = unittest.makeSuite(TaskerTest)
	unittest.TextTestRunner(verbosity=2).run(suite)
