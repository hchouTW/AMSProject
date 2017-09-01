#!/usr/bin/env ruby
require "fileutils"
require "/afs/cern.ch/user/h/hchou/public/BSUB/submit/src/select.rb"
require "/afs/cern.ch/user/h/hchou/public/BSUB/submit/src/parseconfig.rb"

CASTORPrefix = 'xroot://castorpublic.cern.ch/'
CASTORfixOS  = '-OSstagerHost=castorpublic -OSsvcClass=amsuser'
CASTORfixOD  = '-ODstagerHost=castorpublic -ODsvcClass=amsuser'

if (ARGV.length != 0) && (ARGV.length != 2)
	puts "Input arguments is failure. EXITING ..."
	exit
end

if (ARGV.length == 2)
	SUBMODE = "CONFIG"
	ISUI = false
	if (ARGV[0].eql? "-c")
		CONFIG_FILE = ARGV[1]
		if !(File.exist? "#{CONFIG_FILE}")
			puts "#{CONFIG_FILE} isn't exist. EXITING ..."
			exit
		else
			CONFIG = ParseConfig.new(CONFIG_FILE)
		end
	else
		puts "Input arguments should be \'-c\' EXITING ..."
		exit
	end
else
	SUBMODE = "UI"
	ISUI = true
end

SUPER_USER = "hchou"
USER = ENV['USER']
FL = USER.chars.first

HOSTNAME = ENV['HOSTNAME']
IS_LXPLUS = (HOSTNAME.include? "lxplus")

USER_DIR = "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user"

ServiceDir = "#{USER_DIR}/#{USER}/core"
if !(File.exist? ServiceDir)
	puts "#{ServiceDir} isn't exist. EXITING ..."
	exit
end

AMSEnvDir = "/afs/cern.ch/user/h/hchou/public/BSUB/submit/env"
AMSEnv = "setup_amsenv_root6.sh"
if !(File.exist? "#{AMSEnvDir}/#{AMSEnv}")
	puts "#{AMSEnvDir}/#{AMSEnv} isn't exist. EXITING ..."
	exit
end

#LogDir = "/afs/cern.ch/work/#{FL}/#{USER}/" + ((USER.eql? SUPER_USER) ? "" : "public/") + "BSUB_LOG"
LogDir = "/afs/cern.ch/work/#{FL}/#{USER}/BSUB_LOG"

if !(File.exist? LogDir) then Dir.mkdir(LogDir) end

CorePerPCNode = 3
localHost = Array.new

EOScomd = "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select"

# SUBMIT TIME
BsubTime = Time.utc(*(Time.now))
BsubTimeStr = BsubTime.strftime("DATE%Y%m%dCLOCK%H%M%S")

puts "**********************************************************************"
puts "****************************** INPUT PARAMETERS **********************"
QueueModeList = ("localhost" + ((!IS_LXPLUS) ? "" : " ams1nd 1nw 2nw 2nw4cores")).split(" ")
QueueMode = selectShow(SUBMODE, "QUEUE MODE", ((ISUI) ? "" : CONFIG['QUEUE']['MODE']), QueueModeList);

StorageMode = selectShow(SUBMODE, "STORAGE MODE", ((ISUI) ? "" : CONFIG['STORAGE']['MODE']), ((QueueMode.eql? "localhost")? %w(CASTOR EOS NAS AFS) : %w(CASTOR EOS AFS)))

case QueueMode
	when "localhost"
    LocalHost_ListFile = "/afs/cern.ch/user/h/hchou/public/BSUB/submit/src/list.localhost"
		if (File.exist? LocalHost_ListFile)
			HostListMode = (USER.eql? SUPER_USER) ? selectShow(SUBMODE, "HOST LIST MODE", ((ISUI) ? "" : CONFIG["QUEUE"]["HOSTLIST"]), %w(ORIGINAL REBUILD)) : "ORIGINAL"
		else
			HostListMode = "REBUILD"
			if (!(USER.eql? SUPER_USER)); then
	      puts "LocalHist (#{LocalHost_ListFile}) isn't exist. EXITING ..."
			  exit
			end
		end

		case HostListMode
			when "ORIGINAL"
				File.open(LocalHost_ListFile, "r").each(sep="\n") do |line|
					@host = line.chomp
					if (StorageMode.eql? "NAS" and ! @host.include? "pctwams") then next end
					#@CheckComd = "ssh -o StrictHostKeyChecking=no -o BatchMode=yes -o ConnectTimeout=1 #{@host} echo $? &> /dev/null"
					#@result = system(@CheckComd)
					#if ! @result then next end
			  	localHost.push @host
				end
			when "REBUILD"
				@localHost_TW = %w(pctwams02 pctwams05 pctwams06 pctwams07 pctwams10 pctwams11 pctwams12)
				#@localHost_VM = %w(hchou01 hchou02 hchou03 hchou04 hchou05 huli02 zlivm7 zlivm9 zlivm10 zlivm11 zqu01 zqu02 zqu03 zqu04 zqu05 dliuvm01 dliuvm02 dliuvm03 dliuvm04 dliuvm05) # ams.cern.ch
				@localHost_LP = Array.new
				if (! StorageMode.eql? "NAS") then 
					for @index in 000..100
						@hostname = "lxplus%03i" % [@index]
						@checkExist = "host %s.cern.ch &> /dev/null" % [@hostname]
						@result = system(@checkExist)
						if ! @result then next end
						@localHost_LP.push @hostname
					end
					#for @index in 000..100
					#	@hostname = "lxplus%04i" % [@index]
					#	@checkExist = "host %s.cern.ch &> /dev/null" % [@hostname]
					#	@result = system(@checkExist)
					#	if ! @result then next end
					#	@localHost_LP.push @hostname
					#end
				end

				@localHost_CAND = Array.new
				if (StorageMode.eql? "NAS") then 
					@localHost_CAND.push @localHost_TW 
				end
				#@localHost_CAND.push @localHost_VM
				@localHost_CAND.push @localHost_LP
				@localHost_CAND.flatten!

				@localHost_CAND.each { |host|
					@CheckComd = "ssh -o StrictHostKeyChecking=no -o BatchMode=yes -o ConnectTimeout=1 #{host} echo $? &> /dev/null"
					@result = system(@CheckComd)
					if ! @result then next end
					localHost.push host
				}

				File.open(LocalHost_ListFile, "w") { |file|
					file.write localHost.join("\n")
					file.write "\n"
				}
		end
		localHost.shuffle!

		puts "LOCAL HOST (#{localHost.length}) : "
		print localHost.join(' ')
		puts "", ""
	
	when "ams1nd", "1nw", "2nw", "2nw4cores"
		puts "LXPLUS BASH : "
		system("bqueues --q #{QueueMode}")
		puts ""
end

RunMode = selectShow(SUBMODE, "RUN MODE", ((ISUI) ? "" : CONFIG['PROJECT']['RUNMODE']), %w(RUN RERUN))

ProjectMode = selectShow(SUBMODE, "PROJECT MODE", ((ISUI) ? "" : CONFIG['PROJECT']['MODE']), %w(PRODUCTION ANALYSIS))

if (ProjectMode.eql? "ANALYSIS")
  @projDir = ServiceDir + "/analysis"
  if !(File.exist? @projDir)
	  puts "#{@projDir} isn't exist. EXITING ..."
	  exit
  end
  SubProjectMode = selectShow(SUBMODE, "SUB PROJECT MODE", ((ISUI) ? "" : CONFIG['PROJECT']['SUBPROJ']), Dir.entries("#{@projDir}").delete_if{|x| x=="." || x==".."})
end

JobCoreDir = ServiceDir + "/" + ((ProjectMode.eql? "PRODUCTION") ? "production" : "analysis/#{SubProjectMode}")
if !(File.exist? JobCoreDir)
	puts "#{JobCoreDir} isn't exist. EXITING ..."
	exit
end

Version = selectShow(SUBMODE, "VERSION", ((ISUI) ? "" : CONFIG['PROJECT']['VERSION']), Dir.entries("#{JobCoreDir}").delete_if{|x| x=="." || x==".."})

print "VERSION TITLE : "
if (Version.eql? "vdev")
	VersionTitle = BsubTimeStr
else
	if (SUBMODE.eql? "UI")
  	@title = $stdin.gets.chomp
		if ((@tile.eql? "TIME") || (@title.eql? ""))
			VersionTitle = BsubTimeStr
		else
			VersionTitle = @title
		end
	elsif (SUBMODE.eql? "CONFIG")
		@projTitle = CONFIG['PROJECT']['TITLE']
		VersionTitle = (((@projTitle.eql? "TIME") || (@projTitle.eql? "")) ? BsubTimeStr : CONFIG['PROJECT']['TITLE'])
	else
		VersionTitle = BsubTimeStr
	end
end
print "#{VersionTitle}"
puts "", ""

VersionStream = Version + "__" + VersionTitle

JobExeDir = JobCoreDir + "/" + Version
if !(File.exist? JobExeDir)
	puts "#{JobExeDir} isn't exist. EXITING ..."
	exit
else
	puts "JOB EXE FOLDER #{JobExeDir}"
end

#JobExe = selectShow(SUBMODE, "JOB EXE", ((ISUI) ? "" : CONFIG['PROJECT']['JOBEXE']), Dir.entries("#{JobExeDir}").delete_if{|x| x=="." || x==".." || !((x.include? "YiProd") || (x.include? "YiAna"))})
JobExe = selectShow(SUBMODE, "JOB EXE", ((ISUI) ? "" : CONFIG['PROJECT']['JOBEXE']), Dir.entries("#{JobExeDir}").delete_if{ |x| x=="." || x==".." || (x.include? "~") || (x.include? ".conf") || (x.include? ".sh") || (x.include? ".h") || (x.include? ".tcc") || (x.include? ".C") || (x.include? ".so") || (x.include? ".job") })

JobExeProjName = ((((JobExe.split(".")).size()).eql? 0) ? "" : (JobExe.split(".")).at(1))

EventMode = selectShow(SUBMODE, "EVENT MODE", ((ISUI) ? "" : CONFIG['PROJECT']['EVENTMODE']), %w(ISS BT MC))

StreamDir = "#{USER_DIR}/#{USER}/dataset/" + ((ProjectMode.eql? "PRODUCTION") ? "production" : "analysis")

if !(File.exist? StreamDir)
	puts "#{StreamDir} isn't exist. EXITING ..."
	exit
else
	puts "STREAM FOLDER #{StreamDir}"
end
Stream = selectShow(SUBMODE, "STREAM", ((ISUI) ? "" : CONFIG['PROJECT']['STREAM']), Dir.entries("#{StreamDir}").reject { |x| (x.include? "~") || !(x.include? EventMode.downcase) } )
TotalFile = File.open("#{StreamDir+"/"+Stream}", "r").readlines.size

ExeRegionOption = selectShow(SUBMODE, "EXE REGION OPTION (FILE_#{TotalFile})", ((ISUI) ? "" : CONFIG['QUEUE']['REGION']), %w(WHOLE PART))
FilePerExe = inputValue(SUBMODE, "FILE PER EXE", ((ISUI) ? "" : CONFIG['QUEUE']['FILEPEREXE']), [1, TotalFile])
puts ""
case ExeRegionOption
	when "WHOLE"
		ExeStartID = 0
		ExeEndID = (TotalFile / FilePerExe) + ((TotalFile % FilePerExe == 0) ? 0 : 1) - 1
	when "PART"
		@regSExe = 0
		@regEExe = ((TotalFile / FilePerExe) + ((TotalFile % FilePerExe == 0) ? 0 : 1) - 1)
		puts "EXE REGION (#{@regSExe} ... #{@regEExe})"
		ExeStartID = inputValue(SUBMODE, "EXE START ID", ((ISUI) ? "" : CONFIG['QUEUE']['EXESATID']), [@regSExe, @regEExe])
		@regSExe = ExeStartID
		ExeEndID   = inputValue(SUBMODE, "EXE END ID",   ((ISUI) ? "" : CONFIG['QUEUE']['EXEENDID']), [@regSExe, @regEExe])
		puts ""
end
TotalExe = ExeEndID - ExeStartID + 1
puts "EXE REGION (#{TotalExe}) (#{ExeStartID} ... #{ExeEndID})", ""

case QueueMode
	when "ams1nd", "1nw", "2nw", "2nw4cores"
		ExePerJob = inputValue(SUBMODE, "EXE PER JOB", ((ISUI) ? "" : CONFIG['QUEUE']['EXEPERJOB']), [1, TotalExe]);
	when "localhost"
		@TotalPCCore = (localHost.length * CorePerPCNode)
		ExePerJob = (TotalExe / @TotalPCCore) + ((TotalExe % @TotalPCCore == 0) ? 0 : 1)
		puts "EXE PER JOB : #{ExePerJob}"
end
JobStartID = ExeStartID / ExePerJob
JobEndID = ExeEndID / ExePerJob
TotalJob = JobEndID - JobStartID + 1
puts "JOB REGION (#{TotalJob}) (#{JobStartID} ... #{JobEndID})", ""


puts "****************************** OUTPUT PARAMETERS *********************"
puts "SUBMIT_TIME        #{BsubTimeStr}"
puts "QUEUE_MODE         #{QueueMode}"
puts "STORAGE_MODE       #{StorageMode}"
puts "PROJECT_MODE       #{ProjectMode}"
if (ProjectMode.eql? "ANALYSIS"); then
  puts "SUBPROJECT_MODE    #{SubProjectMode}"
end
puts "RUN_MODE           #{RunMode}"
puts "EVENT_MODE         #{EventMode}"
puts ""
puts "VERSION            #{Version}"
puts "VERSION_TITLE      #{VersionTitle}"
puts "VERSION_STREAM     #{VersionStream}"
puts "JOB_EXE            #{JobExe}"
puts "STREAM             #{Stream}"
puts ""
puts "FILE_PER_EXE       #{FilePerExe}"
puts "EXE_PER_JOB        #{ExePerJob}"
puts ""
puts "TOTAL_RUN_EXE      #{TotalExe}"
puts "EXE_START_ID       #{ExeStartID}"
puts "EXE_END_ID         #{ExeEndID}"
puts ""
puts "TOTAL_RUN_JOB      #{TotalJob}"
puts "JOB_START_ID       #{JobStartID}"
puts "JOB_END_ID         #{JobEndID}"
puts "****************************** CONFIRM *******************************"
if (SUBMODE.eql? "UI") || ((SUBMODE.eql? "CONFIG") && (CONFIG['CONFIRM']['MODE'].eql? "NONE"))
	loop do
		print "CONFIRM (YES | NO) : "
		case $stdin.gets.chomp
			when "YES", "yes", "Y", "y"
				puts "START SUBMIT JOB", ""
				break
			when "NO", "no", "N", "n"
				puts "STOP SUBMIT JOB.", ""
				exit
				break
		end
	end
elsif SUBMODE.eql? "CONFIG"
	@confirm = selectShow(SUBMODE, "CONFIRM (YES | NO)", CONFIG['CONFIRM']['MODE'], %w(YES NO));
	case @confirm
		when "YES", "yes", "Y", "y"
			puts "START SUBMIT JOB", ""
		when "NO", "no", "N", "n"
			puts "STOP SUBMIT JOB.", ""
			exit
		else
			puts "Can't find submode. EXITING ...", ""
			exit
	end
else
	puts "Can't find submode. EXITING ...", ""
	exit
end
puts "**********************************************************************"
puts ""
puts "**********************************************************************"
puts "****************************** SUBMIT ********************************"
puts "**********************************************************************"
puts ""


# LOG FOLDER
ProjectLogDir = LogDir + "/" + ProjectMode
if !(File.exist? ProjectLogDir) then Dir.mkdir(ProjectLogDir) end

if (ProjectMode.eql? "PRODUCTION"); then
  StreamLogDir = ProjectLogDir + "/" + Stream
  if !(File.exist? StreamLogDir) then Dir.mkdir(StreamLogDir) end
  VersionStreamLogDir = StreamLogDir + "/" + VersionStream
  if (RunMode.eql? "RUN") && (File.exist? VersionStreamLogDir) then FileUtils.rm_rf VersionStreamLogDir end
else
  SubProjectLogDir = ProjectLogDir + "/" + SubProjectMode
  if !(File.exist? SubProjectLogDir) then Dir.mkdir(SubProjectLogDir) end
  SubProjJobExeLogDir = SubProjectLogDir + ((JobExeProjName.eql? "") ? "" : "/#{JobExeProjName}")
  if !(File.exist? SubProjJobExeLogDir) then Dir.mkdir(SubProjJobExeLogDir) end
  VersionStreamLogDir = SubProjJobExeLogDir + "/" + VersionStream
  if (RunMode.eql? "RUN") && (File.exist? VersionStreamLogDir) then FileUtils.rm_rf VersionStreamLogDir end
end

if !(File.exist? VersionStreamLogDir) then Dir.mkdir(VersionStreamLogDir) end
ExeIDLogDir = VersionStreamLogDir + "/EXEID"
if !(File.exist? ExeIDLogDir) then Dir.mkdir(ExeIDLogDir) end
JobIDLogDir = VersionStreamLogDir + "/JOBID"
if !(File.exist? JobIDLogDir) then Dir.mkdir(JobIDLogDir) end
if QueueMode.eql? "localhost"
	BashJobLogDir = VersionStreamLogDir + "/BASHJOB"
	if !(File.exist? BashJobLogDir) then Dir.mkdir(BashJobLogDir) end
end

if (RunMode.eql? "RUN")
	FileUtils.cp((StreamDir + "/" + Stream), VersionStreamLogDir)
	FileUtils.cp((JobExeDir + "/" + JobExe), VersionStreamLogDir)
	FileUtils.cp((AMSEnvDir + "/" + AMSEnv), VersionStreamLogDir)
end

ParamFileName = "PARAMETERS"
ParamFilePath = "#{VersionStreamLogDir}/#{ParamFileName}"
ParamFile = File.open("#{ParamFilePath}", "w")
ParamFile << """
QUEUE_MODE         #{QueueMode}
STORAGE_MODE       #{StorageMode}
PROJECT_MODE       #{ProjectMode}
"""
if (ProjectMode.eql? "ANALYSIS"); then
ParamFile << """
SUBPROJECT_MODE    #{SubProjectMode}
"""
end
ParamFile << """
RUN_MODE           #{RunMode}
EVENT_MODE         #{EventMode}

VERSION            #{Version}
VERSION_TITLE      #{VersionTitle}
VERSION_STREAM     #{VersionStream}
JOB_EXE            #{JobExe}
STREAM             #{Stream}

FILE_PER_EXE       #{FilePerExe}
EXE_PER_JOB        #{ExePerJob}

TOTAL_RUN_EXE      #{TotalExe}
EXE_START_ID       #{ExeStartID}
EXE_END_ID         #{ExeEndID}

TOTAL_RUN_JOB      #{TotalJob}
JOB_START_ID       #{JobStartID}
JOB_END_ID         #{JobEndID}

"""
ParamFile.close

# OUTPUT FOLDER
case StorageMode
	when "CASTOR"
		StorageServiceDir = "/castor/cern.ch/user/#{FL}/#{USER}/BSUB_#{ProjectMode}"
		%x{nsmkdir -p #{StorageServiceDir}}
    if (ProjectMode.eql? "PRODUCTION"); then
		  StreamSrgDir = StorageServiceDir + "/" + Stream
		  %x{nsmkdir -p #{StreamSrgDir}}
		  VersionStreamSrgDir = StreamSrgDir + "/" + VersionStream
		  %x{nsmkdir -p #{VersionStreamSrgDir}}
    else
		  SubProjectSrgDir = StorageServiceDir + "/" + SubProjectMode
		  %x{nsmkdir -p #{SubProjectSrgDir}}
      SubProjJobExeSrgDir = SubProjectSrgDir + ((JobExeProjName.eql? "") ? "" : "/#{JobExeProjName}")
		  %x{nsmkdir -p #{SubProjJobExeSrgDir}}
      VersionStreamSrgDir = SubProjJobExeSrgDir + "/" + VersionStream
		  %x{nsmkdir -p #{VersionStreamSrgDir}}
    end
	when "EOS"
		StorageServiceDir = "/eos/ams/user/#{FL}/#{USER}/BSUB_#{ProjectMode}"
		%x{#{EOScomd} mkdir -p #{StorageServiceDir}}
    if (ProjectMode.eql? "PRODUCTION"); then
		  StreamSrgDir = StorageServiceDir + "/" + Stream
		  %x{#{EOScomd} mkdir -p #{StreamSrgDir}}
		  VersionStreamSrgDir = StreamSrgDir + "/" + VersionStream
		  %x{#{EOScomd} mkdir -p #{VersionStreamSrgDir}}
    else
      SubProjectSrgDir = StorageServiceDir + "/" + SubProjectMode
		  %x{#{EOScomd} mkdir -p #{SubProjectSrgDir}}
      SubProjJobExeSrgDir = SubProjectSrgDir + ((JobExeProjName.eql? "") ? "" : "/#{JobExeProjName}")
		  %x{#{EOScomd} mkdir -p #{SubProjJobExeSrgDir}}
      VersionStreamSrgDir = SubProjJobExeSrgDir + "/" + VersionStream
		  %x{#{EOScomd} mkdir -p #{VersionStreamSrgDir}}
    end
	when "NAS"
		NASName = "nas06"
		StorageServiceDir = "/#{NASName}/#{USER}/BSUB_#{ProjectMode}"
		if !(File.exist? "/#{NASName}")
			puts "Cannot connect to STORAGE #{NASName}. EXITING ..."
			exit
		else
			if !(File.exist? "/#{NASName}/#{USER}") then Dir.mkdir("/#{NASName}/#{USER}") end
		end
		if !(File.exist? StorageServiceDir) then Dir.mkdir(StorageServiceDir) end
    if (ProjectMode.eql? "PRODUCTION"); then
		  StreamSrgDir = StorageServiceDir + "/" + Stream
		  if !(File.exist? StreamSrgDir) then Dir.mkdir(StreamSrgDir) end
		  VersionStreamSrgDir = StreamSrgDir + "/" + VersionStream
		  if !(File.exist? VersionStreamSrgDir) then Dir.mkdir(VersionStreamSrgDir) end
    else
		  SubProjectSrgDir = StorageServiceDir + "/" + SubProjectMode
		  if !(File.exist? SubProjectSrgDir) then Dir.mkdir(SubProjectSrgDir) end
      SubProjJobExeSrgDir = SubProjectSrgDir + ((JobExeProjName.eql? "") ? "" : "/#{JobExeProjName}")
      if !(File.exist? SubProjJobExeSrgDir) then Dir.mkdir(SubProjJobExeSrgDir) end
      VersionStreamSrgDir = SubProjJobExeSrgDir + "/" + VersionStream
		  if !(File.exist? VersionStreamSrgDir) then Dir.mkdir(VersionStreamSrgDir) end
    end
	when "AFS"
		#StorageServiceDir = "/afs/cern.ch/work/#{FL}/#{USER}/BSUB_#{ProjectMode}"
		StorageServiceDir = "/afs/cern.ch/work/#{FL}/#{USER}/" + ((USER.eql? SUPER_USER) ? "" : "public/") + "BSUB_#{ProjectMode}"
		if !(File.exist? StorageServiceDir) then Dir.mkdir(StorageServiceDir) end
    if (ProjectMode.eql? "PRODUCTION"); then
		  StreamSrgDir = StorageServiceDir + "/" + Stream
		  if !(File.exist? StreamSrgDir) then Dir.mkdir(StreamSrgDir) end
		  VersionStreamSrgDir = StreamSrgDir + "/" + VersionStream
		  if !(File.exist? VersionStreamSrgDir) then Dir.mkdir(VersionStreamSrgDir) end
    else
		  SubProjectSrgDir = StorageServiceDir + "/" + SubProjectMode
		  if !(File.exist? SubProjectSrgDir) then Dir.mkdir(SubProjectSrgDir) end
      SubProjJobExeSrgDir = SubProjectSrgDir + ((JobExeProjName.eql? "") ? "" : "/#{JobExeProjName}")
      if !(File.exist? SubProjJobExeSrgDir) then Dir.mkdir(SubProjJobExeSrgDir) end
      VersionStreamSrgDir = SubProjJobExeSrgDir + "/" + VersionStream
		  if !(File.exist? VersionStreamSrgDir) then Dir.mkdir(VersionStreamSrgDir) end
    end
end

TmpVersionStreamSrgDir= "/tmp/#{USER}/" + VersionStreamSrgDir[VersionStreamSrgDir=~/BSUB_/, VersionStreamSrgDir.length-1].gsub(/\//, "_")

if RunMode.eql? "RUN"
	case StorageMode
		when "CASTOR"
			%x{nsrm -f #{VersionStreamSrgDir}/*}
		when "EOS"
			%x{#{EOScomd} rm -f #{VersionStreamSrgDir}/*}
		when "NAS"
			%x{/bin/rm -f #{VersionStreamSrgDir}/*}
		when "AFS"
			%x{/bin/rm -f #{VersionStreamSrgDir}/*}
	end
end
	
# Copy PARAMETERS To StorageDir
case StorageMode
	when "CASTOR"
		#%x{rfcp #{ParamFilePath} #{VersionStreamSrgDir}/#{ParamFileName}}
		%x{xrdcp -f #{ParamFilePath} #{CASTORPrefix}/#{VersionStreamSrgDir}/#{ParamFileName} #{CASTORfixOD}}
	when "EOS"
		%x{#{EOScomd} cp #{ParamFilePath} #{VersionStreamSrgDir}/#{ParamFileName}}
	when "NAS"
		%x{/bin/cp #{ParamFilePath} #{VersionStreamSrgDir}/#{ParamFileName}}
	when "AFS"
		%x{/bin/cp #{ParamFilePath} #{VersionStreamSrgDir}/#{ParamFileName}}
end

##############################################################################
# JobScript
JobScript = "JobScript.sh"
JobScriptFile = File.open("#{VersionStreamLogDir}/#{JobScript}", "w")
JobScriptFile << "#!bin/bash\n"
JobScriptFile << "#shopt -s -o nounset\n"
JobScriptFile << "CASTORPrefix='xroot://castorpublic.cern.ch/'\n"
JobScriptFile << "CASTORfixOS='-OSstagerHost=castorpublic -OSsvcClass=amsuser'\n"
JobScriptFile << "CASTORfixOD='-ODstagerHost=castorpublic -ODsvcClass=amsuser'\n"
JobScriptFile << "EOScomd='/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'\n"
JobScriptFile << """
echo -e \"====  Start Run Time : '$(date)'  ====\"
echo -e \"====  Local Host : "'$HOSTNAME'"  ====\"
echo -e \"====  Redhat-release  ====\"
cat /etc/redhat-release

echo -e \"*************************\"
echo -e \"****  START RUNNING  ****\"
echo -e \"*************************\"

declare -i JobID
declare -i ExeStartID
declare -i ExeEndID

JobID=$1
ExeStartID=$2
ExeEndID=$3

"""

case QueueMode
	when "ams1nd", "1nw", "2nw", "2nw4cores"
JobScriptFile << """
TmpVersionStreamSrgDir=$(pwd)
"""
	when "localhost"
JobScriptFile << """
TmpVersionStreamSrgDir=#{TmpVersionStreamSrgDir}_JOBID${JobID}
mkdir -p ${TmpVersionStreamSrgDir}

"""
end

JobScriptFile << """
VersionStreamSrgDir=#{VersionStreamSrgDir}
VersionStreamLogDir=#{VersionStreamLogDir}
ExeIDLogDir=#{ExeIDLogDir}
Stream=#{Stream}
"""

case QueueMode
	when "ams1nd", "1nw", "2nw", "2nw4cores"
JobScriptFile << """
cp ${VersionStreamLogDir}/#{Stream} .
cp ${VersionStreamLogDir}/#{JobExe} .
cp ${VersionStreamLogDir}/#{AMSEnv} .
ls -alth
echo -e \"source $(pwd)/#{AMSEnv} icc64 root604\"
source $(pwd)/#{AMSEnv} icc64 root604
echo $LD_LIBRARY_PATH
ldd #{JobExe}
"""
	when "localhost"
JobScriptFile << """
ls -alth ${VersionStreamLogDir}
echo -e \"source ${VersionStreamLogDir}/#{AMSEnv}\"
source ${VersionStreamLogDir}/#{AMSEnv} icc64 root604
echo $LD_LIBRARY_PATH
ldd ${VersionStreamLogDir}/#{JobExe}
"""
end

JobScriptFile << """

for (( exeID=${ExeStartID}; exeID<=${ExeEndID}; exeID++ ))
do
	logID=$(printf \"JOB_%07i_EXE_%07i\" $JobID $exeID)
	logID_File=${ExeIDLogDir}/${logID}
"""
if RunMode.eql? "RERUN"
JobScriptFile << """
	# check old run is finished.
	if [ -f ${logID_File} ]; then
		isFind=`grep -L \"****  FINISH RUNNING  ****\" ${logID_File}`
		if [[ ${isFind} != ${logID_File} ]]; then
			continue
		else
			/bin/rm ${logID_File}
		fi
	fi
"""
end

JobScriptFile << """
	echo -e \"****  RUNNING  ****\"
	echo -e \"****  RUNNING  ****\" | tee ${logID_File}
"""

case QueueMode
	when "ams1nd", "1nw", "2nw", "2nw4cores"
JobScriptFile << """
	echo -e \"./#{JobExe} #{EventMode} ${Stream} ${exeID} #{FilePerExe} ${TmpVersionStreamSrgDir}\"
	./#{JobExe} #{EventMode} ${Stream} ${exeID} #{FilePerExe} ${TmpVersionStreamSrgDir} 2>&1 | tee -a ${logID_File}
"""
	when "localhost"
JobScriptFile << """
	echo -e \"${VersionStreamLogDir}/#{JobExe} #{EventMode} ${VersionStreamLogDir}/${Stream} ${exeID} #{FilePerExe} ${TmpVersionStreamSrgDir}\"
	${VersionStreamLogDir}/#{JobExe} #{EventMode} ${VersionStreamLogDir}/${Stream} ${exeID} #{FilePerExe} ${TmpVersionStreamSrgDir} 2>&1 | tee -a ${logID_File}
"""
end

JobScriptFile << """
	ls -alth ${TmpVersionStreamSrgDir}
	fileName=$(printf \"%07i.root\" $exeID)
	rootFile=`ls ${TmpVersionStreamSrgDir} | grep ${fileName}`
"""

case StorageMode
	when "CASTOR"
		#JobScriptFile << "	rfcp ${TmpVersionStreamSrgDir}/$rootFile ${VersionStreamSrgDir}/$rootFile"
		JobScriptFile << "	xrdcp -f ${TmpVersionStreamSrgDir}/$rootFile ${CASTORPrefix}/${VersionStreamSrgDir}/$rootFile ${CASTORfixOD}"
	when "EOS"
		JobScriptFile << "	$EOScomd cp ${TmpVersionStreamSrgDir}/$rootFile ${VersionStreamSrgDir}/$rootFile"
	when "NAS"
		JobScriptFile << "	cp -f ${TmpVersionStreamSrgDir}/$rootFile ${VersionStreamSrgDir}/$rootFile"
	when "AFS"
		JobScriptFile << "	cp -f ${TmpVersionStreamSrgDir}/$rootFile ${VersionStreamSrgDir}/$rootFile"
end

JobScriptFile << """
	/bin/rm ${TmpVersionStreamSrgDir}/$rootFile

	echo -e \"****  FINISH RUNNING  ****\" | tee -a ${logID_File}
	echo -e \"****  FINISH RUNNING  ****\"
	echo -e \"\\n\"
done

"""

if QueueMode.eql? "localhost"
JobScriptFile << """
/bin/rm -rf ${TmpVersionStreamSrgDir}

"""
end

JobScriptFile << """
echo -e \"**************************\"
echo -e \"****  FINISH RUNNING  ****\"
echo -e \"**************************\"
"""
JobScriptFile.close
##############################################################################

## Submit Job
for @jobID in JobStartID..JobEndID
	@exeSID = @jobID * ExePerJob
	if @exeSID < ExeStartID then @exeSID = ExeStartID end
	@exeEID = (@jobID + 1) * ExePerJob - 1
	if @exeEID > ExeEndID then @exeEID = ExeEndID end
	@CountJobs = @jobID - JobStartID + 1

	@BsubName = ((ProjectMode.eql? "PRODUCTION") ? "#{ProjectMode}__#{VersionStream}__#{Stream}__JOBID%07i" : "#{ProjectMode}__#{SubProjectMode}__#{VersionStream}__JOBID%07i") % [@jobID]
	@LogFile = "#{JobIDLogDir}/log.#{@BsubName}"

	# check old run is finished.
	if RunMode.eql? "RERUN"
		@isOK = true
		for @exeID in @exeSID..@exeEID
			@logID = "JOB_%07i_EXE_%07i" % [@jobID, @exeID]
			@logID_File = "#{ExeIDLogDir}/#{@logID}"
			if (File.exist? @logID_File)
				@isOK = File.read(@logID_File).include?("****  FINISH RUNNING  ****")
				if @isOK then next end
			else
				@isOK = false
			end
			if !@isOK then break end
		end
		if @isOK then next end
	end

	case QueueMode
		when "ams1nd", "1nw", "2nw", "2nw4cores"
			@BSUB_JobComd = "bsub -q #{QueueMode} -J #{@BsubName} -oo #{@LogFile} sh #{VersionStreamLogDir}/#{JobScript} #{@jobID} #{@exeSID} #{@exeEID}"
			system("#{@BSUB_JobComd}")
		when "localhost"
			@PCNode = (@CountJobs -1) % (localHost.length)
			@PCName = localHost[@PCNode]
			@SSH_JobComd = "ssh -o StrictHostKeyChecking=no -o BatchMode=yes -o ConnectTimeout=1 #{@PCName} sh #{VersionStreamLogDir}/#{JobScript} #{@jobID} #{@exeSID} #{@exeEID} &> #{@LogFile} &"

			@bashJob = "#{BashJobLogDir}/bashJob%07i.sh" % @jobID
			@bashJobFile = File.open("#{@bashJob}", "w")
			@bashJobFile << "#!bin/bash\n"
			@bashJobFile << "#{@SSH_JobComd}\n"
			@bashJobFile.close

			puts "#{@BsubName}  #{@PCName}"
			system("sh #{@bashJob}")
	end
end

puts ""
puts "**********************************************************************"
puts "************************** SUBMIT FINISH *****************************"
puts "**********************************************************************"
puts ""
