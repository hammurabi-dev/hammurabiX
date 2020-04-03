#ifndef HAMMURABI_IO_H
#define HAMMURABI_IO_H

#include <cassert>
#include <fstream>
#include <hamdis.h>
#include <hamp.h>
#include <hamtype.h>
#include <iostream>
#include <memory>
#include <omp.h>
#include <stdexcept>
#include <string>

template <typename T> class Hamio {
protected:
  std::string Filename;

public:
  Hamio() = default;
  Hamio(const std::string &filename) { this->Filename = filename; }
  Hamio &operator=(const Hamio<T> &) = delete;
  Hamio &operator=(Hamio<T> &&) = delete;
  Hamio(const Hamio<T> &) = delete;
  Hamio(Hamio<T> &&) = delete;
  virtual ~Hamio() = default;
  // filename
  virtual std::string filename() const { return this->Filename; }
  // set filename
  virtual void filename(const std::string &filename) {
    this->Filename = filename;
  }
  // write to disk (for Hamis)
  virtual void dump(const Hamdis<T> &m) const {
    std::fstream outfile(this->Filename.c_str(),
                         std::ios::out | std::ios::binary);
    if (outfile.is_open()) {
      ham_float tmpfloat;
      ham_uint tmpuint;
      // unable to multithread
      for (ham_uint i = 0; i < m.npix(); ++i) {
        // idx
        tmpuint = m.index(i);
        outfile.write(reinterpret_cast<char *>(&tmpuint), sizeof(ham_uint));
        // theta
        tmpfloat = m.pointing(i).theta();
        outfile.write(reinterpret_cast<char *>(&tmpfloat), sizeof(ham_float));
        // phi
        tmpfloat = m.pointing(i).phi();
        outfile.write(reinterpret_cast<char *>(&tmpfloat), sizeof(ham_float));
        // data
        tmpfloat = m.data(i);
        outfile.write(reinterpret_cast<char *>(&tmpfloat), sizeof(ham_float));
      }
      outfile.close();
    } else {
      std::cerr << this->Filename << std::endl;
      throw std::runtime_error("unable to open");
    }
  }
  // write to disk (for Hampix)
  virtual void dump(const Hampix<T> &m) const {
    std::fstream outfile(this->Filename.c_str(),
                         std::ios::out | std::ios::binary);
    if (outfile.is_open()) {
      ham_float tmpfloat;
      // unable to multithread
      for (ham_uint i = 0; i < m.npix(); ++i) {
        // data
        tmpfloat = m.data(i);
        outfile.write(reinterpret_cast<char *>(&tmpfloat), sizeof(ham_float));
      }
      outfile.close();
    } else {
      std::cerr << this->Filename << std::endl;
      throw std::runtime_error("unable to open");
    }
  }
  // read from disk (for Hamdis)
  virtual void load(Hamdis<T> &m) const {
    std::fstream infile(this->Filename.c_str(),
                        std::ios::in | std::ios::binary);
    if (infile.is_open()) {
      T tmpfloat, tmpfloat2;
      ham_uint tmpuint;
      // unable to multithread
      for (ham_uint i = 0; i < m.npix(); ++i) {
        if (infile.eof())
          throw std::runtime_error("unexpected end of file");
        // idx
        infile.read(reinterpret_cast<char *>(&tmpuint), sizeof(ham_uint));
        m.index(tmpuint, tmpuint);
        // theta
        infile.read(reinterpret_cast<char *>(&tmpfloat), sizeof(T));
        // phi
        infile.read(reinterpret_cast<char *>(&tmpfloat2), sizeof(T));
        m.pointing(tmpuint, tmpfloat, tmpfloat2);
        // data
        infile.read(reinterpret_cast<char *>(&tmpfloat), sizeof(T));
        m.data(tmpuint, tmpfloat);
      }
      infile.close();
    } else {
      std::cerr << this->Filename << std::endl;
      throw std::runtime_error("unable to open");
    }
  }
  // read from disk (for Hampix)
  virtual void load(Hampix<T> &m) const {
    std::fstream infile(this->Filename.c_str(),
                        std::ios::in | std::ios::binary);
    if (infile.is_open()) {
      T tmpfloat;
      // unable to multithread
      for (ham_uint i = 0; i < m.npix(); ++i) {
        if (infile.eof())
          throw std::runtime_error("unexpected end of file");
        // data
        infile.read(reinterpret_cast<char *>(&tmpfloat), sizeof(T));
        m.data(i, tmpfloat);
      }
      infile.close();
    } else {
      std::cerr << this->Filename << std::endl;
      throw std::runtime_error("unable to open");
    }
  }
};

#endif
