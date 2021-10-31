#!/usr/bin/env groovy
import static nextflow.Nextflow.file
import nextflow.Channel

class LikeletUtils {

  // adjust command colors
  static String ANSI_RESET = "\u001B[0m"
  static String ANSI_BLACK = "\u001B[30m"
  static String ANSI_RED = "\u001B[31m"
  static String ANSI_GREEN = "\u001B[32m"
  static String ANSI_YELLOW = "\u001B[33m"
  static String ANSI_BLUE = "\u001B[34m"
  static String ANSI_PURPLE = "\u001B[35m"
  static String ANSI_CYAN = "\u001B[36m"
  static String ANSI_WHITE = "\u001B[37m"

  static def print_red = {  str -> LikeletUtils.ANSI_RED + str + LikeletUtils.ANSI_RESET }
  static def print_black = {  str -> LikeletUtils.ANSI_BLACK + str + LikeletUtils.ANSI_RESET }
  static def print_green = {  str -> LikeletUtils.ANSI_GREEN + str + LikeletUtils.ANSI_RESET }
  static def print_yellow = {  str -> LikeletUtils.ANSI_YELLOW + str + LikeletUtils.ANSI_RESET }
  static def print_blue = {  str -> LikeletUtils.ANSI_BLUE + str + LikeletUtils.ANSI_RESET }
  static def print_cyan = {  str -> LikeletUtils.ANSI_CYAN + str + LikeletUtils.ANSI_RESET }
  static def print_purple = {  str -> LikeletUtils.ANSI_PURPLE + str + LikeletUtils.ANSI_RESET }
  static def print_white = {  str -> LikeletUtils.ANSI_WHITE + str + LikeletUtils.ANSI_RESET }

  // Check if a row has the expected number of item, adjusted from Sarek 
    static def checkNumberOfItem(row, number) {
      if (row.size() != number) exit 1, println("Malformed row in TSV file: ${row}, see --help for more information")
      return true
    } 

  // Return status [0,1]
    // 0 == Normal, 1 == Tumor
    static def returnStatus(it) {
      if (!(it in [0, 1])) exit 1, println("Status is not recognized in TSV file: ${it}, see --help for more information")
      return it
    }

    // Return file if it exists
    static def returnFile(it) {
      if (!file(it).exists()) exit 1, println("Missing file in TSV file: ${it}, see --help for more information")
      return file(it)
    }

  static def sysucc_ascii() {
    print LikeletUtils.print_yellow(" ▄▄▄▄▄▄▄▄▄▄▄  ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ \n")
    print LikeletUtils.print_yellow("▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌\n")
    print LikeletUtils.print_yellow("▐░█▀▀▀▀▀▀▀▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ \n")
    print LikeletUtils.print_yellow("▐░▌          ▐░▌       ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌          ▐░▌          \n")
    print LikeletUtils.print_yellow("▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌▐░▌          ▐░▌          \n")
    print LikeletUtils.print_yellow("▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌          ▐░▌          \n")
    print LikeletUtils.print_yellow(" ▀▀▀▀▀▀▀▀▀█░▌ ▀▀▀▀█░█▀▀▀▀  ▀▀▀▀▀▀▀▀▀█░▌▐░▌       ▐░▌▐░▌          ▐░▌          \n")
    print LikeletUtils.print_yellow("          ▐░▌     ▐░▌               ▐░▌▐░▌       ▐░▌▐░▌          ▐░▌          \n")
    print LikeletUtils.print_yellow(" ▄▄▄▄▄▄▄▄▄█░▌     ▐░▌      ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ \n")
    print LikeletUtils.print_yellow("▐░░░░░░░░░░░▌     ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌\n")
    print LikeletUtils.print_yellow(" ▀▀▀▀▀▀▀▀▀▀▀       ▀       ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ \n")
    }
  // extrct fastq information from tsvFile
  static def extractData(tsvFile) {
    // Channeling the TSV file containing FASTQ.
    // Format is: "subject gender status sample lane fastq1 fastq2"
    def inputData = Channel.from(tsvFile)
    .splitCsv(sep: '\t', skip: 1)
    .map { row ->
      LikeletUtils.checkNumberOfItem(row, 6)
      def idSample  = row[0]
      def fastqFile1 = file(row[1])
      def fastqFile2 = file(row[2])
      def group = row[5]
      def input = true
      def gzip = false
      def readsSingle = false
      def filetype = "fastq"
      if (row[1].endsWith(".gz") == true ){
        gzip = true
      }else if (row[1].endsWith(".bam") == true ){
        filetype = "bam"
      }
      if (row[2].endsWith("false") == true){
        readsSingle = true
        [idSample, [fastqFile1], readsSingle, gzip, input, group, filetype]
      } else {
        [idSample, [fastqFile1, fastqFile2], readsSingle, gzip, input, group, filetype]
      }
    }
    def ipData = Channel.from(tsvFile)
    .splitCsv(sep: '\t', skip: 1)
    .map { row ->
      LikeletUtils.checkNumberOfItem(row, 6)
      def idSample  = row[0]
      def fastqFile1 = file(row[3])
      def fastqFile2 = file(row[4])
      def group = row[5]
      def input = false
      def gzip = false
      def readsSingle = false
      def filetype = "fastq"
      if (row[3].endsWith(".gz") == true){
        gzip = true
      }else if (row[3].endsWith(".bam") == true){
        filetype = "bam"
      }
      if (row[4].endsWith("false") == true){
        readsSingle = true
        [idSample, [fastqFile1], readsSingle, gzip, input, group, filetype]
      } else {
        [idSample, [fastqFile1, fastqFile2], readsSingle, gzip, input, group, filetype]
      }          
    }
    return inputData.mix(ipData)
  }
  static def addstringToalign(String str,int num){
    if(str.length() < num) {
      def numSpace =  num-str.length() 

      numSpace.times{
        str += ' '
      }
    }
    str
  }
}


